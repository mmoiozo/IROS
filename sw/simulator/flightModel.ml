open Stdlib

type meter = float
type meter_s = float
type radian = float
type radian_s = float
type state = {
    start : float;
    mutable t : float;
    mutable x : meter;
    mutable y : meter;
    mutable z : meter;
    mutable psi : radian; (* Trigonometric *)
    mutable phi : radian;
    mutable phi_dot : radian_s;
    mutable delta_a : float;
    mutable thrust : float;
    mutable air_speed : meter_s
  }

let init route = {
  start = Unix.gettimeofday (); t = 0.; x = 0.; y = 0. ; z = 0.;
  psi = route; phi = 0.; phi_dot = 0.;
  delta_a = 0.; thrust = 0.; air_speed = 0.
}

let get_xyz state = (state.x, state.y, state.z)
let get_time state = state.t
let get_phi state = state.phi

let set_air_speed state x = state.air_speed <- x

let drag = 0.45
let c_lp = -10.
let g = 9.81
let weight = 1. *. 1.4

let max_phi = 0.7 (* rad *)
let bound = fun x mi ma -> if x > ma then ma else if x < mi then mi else x



module Make(A:Data.MISSION) = struct
  open Data
(* Minimum complexity *)
(*
   http://controls.ae.gatech.edu/papers/johnson_dasc_01.pdf
   http://controls.ae.gatech.edu/papers/johnson_mst_01.pdf
 *)

  let state_update = fun state (wx, wy) ->
    let now = Unix.gettimeofday () -. state.start in
    let dt = now -. state.t in
    if state.air_speed > 0. then begin
      let phi_dot_dot = state.delta_a +. c_lp *. state.phi_dot /. state.air_speed in
      state.phi_dot <- state.phi_dot +. phi_dot_dot *. dt;
      state.phi <- bound (state.phi +. state.phi_dot *. dt) (-.max_phi) max_phi;
      let psi_dot = -. g /. state.air_speed *. tan state.phi in
      state.psi <- norm_angle (state.psi +. psi_dot *. dt);
      let dx = state.air_speed *. cos state.psi *. dt +. wx *. dt
      and dy = state.air_speed *. sin state.psi *. dt +. wy *. dt in
      state.x <- state.x +.dx ;
      state.y <- state.y +. dy;
      let gamma = (state.thrust -. drag) /. weight in
      let dz = sin gamma *. state.air_speed *. dt in
      state.z <- state.z +. dz
    end;
    state.t <- now


  let servos =
    try
      ExtXml.child A.ac.airframe "servos"
    with
      Not_found ->
	failwith (Printf.sprintf "Child 'servos' expected in '%s'\n" (Xml.to_string A.ac.airframe))
	  
  let misc_section =
    try
      ExtXml.child A.ac.airframe ~select:(fun x -> ExtXml.attrib x "name" = "MISC") "section"
    with
      Not_found ->
	failwith (Printf.sprintf "Child 'section' with 'name=MISC' expected in '%s'\n" (Xml.to_string A.ac.airframe))

  let nominal_airspeed =
    try
      float_of_string (Xml.attrib (ExtXml.child misc_section ~select:(fun x -> ExtXml.attrib x "name" = "NOMINAL_AIRSPEED") "define") "value")
    with
      Not_found ->
	failwith (Printf.sprintf "Child 'define' with 'name=NOMINAL_AIRSPEED' expected in '%s'\n" (Xml.to_string misc_section))

  let get_servo name =
    try
      ExtXml.child servos ~select:(fun x -> ExtXml.attrib x "name" = name) "servo"
    with
      Not_found ->
	failwith (Printf.sprintf "Child 'servo' with name='%s' expected in '%s'\n" name (Xml.to_string servos))

  let us_attrib = fun x a -> int_of_string (ExtXml.attrib x a)

  let gaz = get_servo "GAZ"
  let min_thrust =  us_attrib gaz "min"
  let max_thrust =  us_attrib gaz "max"
      
  type servo_id = int
  type ppm = int

  let no_thrust = int_of_string (ExtXml.attrib gaz "no")

      

  let some_aileron_left = try Some (get_servo "AILERON_LEFT") with _ -> None
  let some_ailevon_left = try Some (get_servo "AILEVON_LEFT") with _ -> None
  let some_ailevon_right = try Some (get_servo "AILEVON_RIGHT") with _ -> None

  let float_attrib = fun x a -> float_of_string (ExtXml.attrib x a)
  let int_attrib = fun x a -> int_of_string (ExtXml.attrib x a)

  let sign = fun x ->
    if float_attrib x "min" < float_attrib x "max" then 1 else -1

  let do_thrust = fun state servo ->
    state.thrust <- (float (servo.(no_thrust) - min_thrust) /. float (max_thrust - min_thrust))

  let do_servos =  
    match some_aileron_left, some_ailevon_left, some_ailevon_right with
      Some aileron_left, None, None ->
	let c_lda = 16. *. 9e-5 in (* phi_dot_dot from aileron *)

	let sign_aileron_left = sign aileron_left
	and n_delta_a = us_attrib aileron_left "neutral"
	and no_aileron_left = int_attrib aileron_left "no" in
	fun state servo ->
	  (** Printf.printf "left=%d\n" (servo.(no_aileron_left) - n_delta_a); flush stdout;  **)
	  state.delta_a <- c_lda *. float (- sign_aileron_left * (servo.(no_aileron_left) - n_delta_a));
	  do_thrust state servo
    | None, Some ailevon_left, Some ailevon_right ->
	let c_lda = 2.5e-4 in (* phi_dot_dot from aileron *)

	let sign_ailevon_left = sign ailevon_left
	and sign_ailevon_right = sign ailevon_right
	and left_neutral = us_attrib ailevon_left "neutral"
	and right_neutral = us_attrib ailevon_right "neutral"
	and left_travel = float (us_attrib ailevon_left "max" - us_attrib ailevon_left "min") /. 1200.
	and right_travel = float (us_attrib ailevon_right "max" - us_attrib ailevon_right "min") /. 1200.
	and no_ailevon_left = int_attrib ailevon_left "no"
	and no_ailevon_right = int_attrib ailevon_right "no" in
	fun state servo ->
	  do_thrust state servo;
	  let sum = (float (servo.(no_ailevon_left) - left_neutral) /. left_travel +.
		       float (servo.(no_ailevon_right) - right_neutral) /.  right_travel) /. 2. in
(*	Printf.printf "%d %f\n" (servo no_ailevon_left - left_neutral) sum; flush stdout; *)
	  state.delta_a <- c_lda *. (-. sum)
    | _ -> failwith "Aileron or Ailevon left and right PLEASE"

  let nb_servos = 10 (* 4017 *)
end

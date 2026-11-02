class_name VesselView
extends Node3D
## The biphasic vessel in three dimensions, read from a recording: Buhse's 250 mL round-bottom flask, its water and
## its oil at their recorded volumes (to scale, as measured after decanting), the mixed aggregates and the empty
## micelles in the water as points (one point for the same number of objects, the larger stock 300 points), the
## magnetic bar that stirs it; over it, the loop with its flows, the numbers of the frame and the curves.

const RADIUS := 0.965  # the liquid's sphere, inside a glass of radius 1
const MOST_POINTS := 300

var recording: Recording
var _liquid: ShaderMaterial
var _surface: MeshInstance3D
var _boundary: MeshInstance3D
var _aggregates: MultiMeshInstance3D
var _micelles: MultiMeshInstance3D
var _bar: MeshInstance3D
var _loop: VesselLoop
var _chart: VesselChart
var _numbers: VBoxContainer
var _water_top: float = 0.0
var _seeds: Array[Vector3] = []
var _counts := Vector2i.ZERO


func _init(source: Recording) -> void:
	recording = source


func _ready() -> void:
	_build_stage()
	_build_flask()
	_build_overlay()


func _build_stage() -> void:
	var environment := Environment.new()
	environment.background_mode = Environment.BG_COLOR
	environment.background_color = Look.NIGHT
	environment.ambient_light_source = Environment.AMBIENT_SOURCE_COLOR
	environment.ambient_light_color = Color(0.55, 0.62, 0.75)
	environment.ambient_light_energy = 0.6
	environment.tonemap_mode = Environment.TONE_MAPPER_FILMIC
	environment.glow_enabled = true
	environment.glow_intensity = 0.4
	var world := WorldEnvironment.new()
	world.environment = environment
	add_child(world)
	var sun := DirectionalLight3D.new()
	sun.rotation_degrees = Vector3(-50, 35, 0)
	sun.light_energy = 1.3
	sun.shadow_enabled = true
	add_child(sun)
	var fill := OmniLight3D.new()
	fill.position = Vector3(-2.2, 0.6, 2.0)
	fill.light_color = Color(1.0, 0.86, 0.66)
	fill.light_energy = 1.2
	fill.omni_range = 6.0
	add_child(fill)
	var bench := MeshInstance3D.new()
	var plane := PlaneMesh.new()
	plane.size = Vector2(12, 12)
	bench.mesh = plane
	bench.position = Vector3(0, -1.02, 0)
	bench.material_override = _matte(Color(0.16, 0.19, 0.25))
	add_child(bench)
	var camera := OrbitCamera.new()
	camera.target = Vector3(0, 0.1, 0)
	camera.distance = 4.2
	add_child(camera)


func _build_flask() -> void:
	var glass := StandardMaterial3D.new()
	glass.albedo_color = Color(0.85, 0.93, 1.0, 0.13)
	glass.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA
	glass.roughness = 0.04
	glass.metallic_specular = 0.9
	glass.rim_enabled = true
	glass.rim = 0.6
	var body := MeshInstance3D.new()
	var sphere := SphereMesh.new()
	sphere.radius = 1.0
	sphere.height = 2.0
	sphere.radial_segments = 96
	sphere.rings = 48
	body.mesh = sphere
	body.material_override = glass
	add_child(body)
	var neck := MeshInstance3D.new()
	var tube := CylinderMesh.new()
	tube.top_radius = 0.22
	tube.bottom_radius = 0.22
	tube.height = 0.9
	tube.cap_top = false
	tube.cap_bottom = false
	neck.mesh = tube
	neck.position = Vector3(0, 1.4, 0)
	neck.material_override = glass
	add_child(neck)
	_liquid = ShaderMaterial.new()
	_liquid.shader = load("res://shaders/liquid.gdshader") as Shader
	var liquid := MeshInstance3D.new()
	var drop := SphereMesh.new()
	drop.radius = RADIUS
	drop.height = 2.0 * RADIUS
	drop.radial_segments = 96
	drop.rings = 48
	liquid.mesh = drop
	liquid.material_override = _liquid
	add_child(liquid)
	_surface = _disc(Color(0.90, 0.69, 0.26, 0.95))
	_boundary = _disc(Color(0.45, 0.62, 0.80, 0.6))
	_aggregates = _points(Look.LIVING, 0.024)
	_micelles = _points(Look.VIOLET, 0.019)
	_bar = MeshInstance3D.new()
	var capsule := CapsuleMesh.new()
	capsule.radius = 0.045
	capsule.height = 0.42
	_bar.mesh = capsule
	_bar.rotation_degrees = Vector3(0, 0, 90)
	_bar.position = Vector3(0, -RADIUS + 0.06, 0)
	_bar.material_override = _matte(Color(0.95, 0.96, 0.98))
	add_child(_bar)


func _build_overlay() -> void:
	var layer := CanvasLayer.new()
	add_child(layer)
	var title := VBoxContainer.new()
	title.position = Vector2(40, 84)
	layer.add_child(title)
	title.add_child(Look.label(recording.title, 30, Color.WHITE, "title"))
	var subtitle := Look.label(recording.subtitle, 15, Look.NIGHT_MUTED)
	title.add_child(subtitle)
	title.add_child(Look.label("Stirred at 800 rpm, 80 °C · volumes to the scale of a 250 mL flask", 13, Look.NIGHT_MUTED))
	var numbers := PanelContainer.new()
	numbers.add_theme_stylebox_override("panel", Look.panel(Color(1, 1, 1, 0.93), 16, 18))
	numbers.position = Vector2(40, 600)
	numbers.custom_minimum_size = Vector2(380, 0)
	layer.add_child(numbers)
	_numbers = VBoxContainer.new()
	_numbers.add_theme_constant_override("separation", 6)
	numbers.add_child(_numbers)
	_loop = VesselLoop.new()
	_loop.position = Vector2(1180, 84)
	layer.add_child(_loop)
	_chart = VesselChart.new()
	_chart.position = Vector2(1180, 404)
	layer.add_child(_chart)


func _matte(colour: Color) -> StandardMaterial3D:
	var material := StandardMaterial3D.new()
	material.albedo_color = colour
	material.roughness = 0.8
	return material


func _disc(colour: Color) -> MeshInstance3D:
	var disc := MeshInstance3D.new()
	var cylinder := CylinderMesh.new()
	cylinder.height = 0.004
	cylinder.radial_segments = 96
	disc.mesh = cylinder
	var material := StandardMaterial3D.new()
	material.albedo_color = colour
	material.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA
	material.roughness = 0.2
	disc.material_override = material
	add_child(disc)
	return disc


func _points(colour: Color, radius: float) -> MultiMeshInstance3D:
	var dot := SphereMesh.new()
	dot.radius = radius
	dot.height = 2.0 * radius
	dot.radial_segments = 10
	dot.rings = 6
	var material := StandardMaterial3D.new()
	material.albedo_color = colour
	material.emission_enabled = true
	material.emission = colour
	material.emission_energy_multiplier = 0.6
	dot.material = material
	var many := MultiMesh.new()
	many.transform_format = MultiMesh.TRANSFORM_3D
	many.mesh = dot
	many.instance_count = MOST_POINTS
	many.visible_instance_count = 0
	var instance := MultiMeshInstance3D.new()
	instance.multimesh = many
	add_child(instance)
	return instance


## The height of liquid in the flask's sphere holding `millilitres` of the 250 mL, in the sphere's coordinates.
func _level(millilitres: float) -> float:
	var fraction := clampf(millilitres / recording.constant("flask_ml", 250.0), 0.0, 1.0)
	var low := 0.0
	var high := 2.0
	for i in 50:
		var h := 0.5 * (low + high)
		if h * h * (3.0 - h) / 4.0 < fraction:
			low = h
		else:
			high = h
	return -RADIUS + 0.5 * (low + high) * RADIUS


func show_frame(index: int) -> void:
	var frame := recording.frame(index)
	if frame == null:
		return
	var water := recording.value(frame, "water_ml")
	var interface := recording.value(frame, "interface_ml")
	var oil := recording.value(frame, "oil_ml")
	_water_top = _level(water)
	var oil_top := _level(water + interface + oil)
	_liquid.set_shader_parameter("water_top", _water_top)
	_liquid.set_shader_parameter("oil_top", oil_top)
	_place_disc(_surface, oil_top, oil_top - _water_top > 0.002)
	_place_disc(_boundary, _water_top, true)
	var mixed := recording.value(frame, "mixed_mol")
	var micelles := recording.value(frame, "micelles_mol")
	var per_point := maxf(maxf(mixed, micelles) / MOST_POINTS, 1e-6)
	_counts = Vector2i(int(minf(mixed / per_point, MOST_POINTS)), int(minf(micelles / per_point, MOST_POINTS)))
	_scatter()
	_loop.show_frame(recording, frame)
	_chart.show_frame(recording, frame)
	_show_numbers(frame)


func _place_disc(disc: MeshInstance3D, height: float, visible_now: bool) -> void:
	var cylinder := disc.mesh as CylinderMesh
	var radius := sqrt(maxf(RADIUS * RADIUS - height * height, 0.0))
	cylinder.top_radius = radius
	cylinder.bottom_radius = radius
	disc.position = Vector3(0, height, 0)
	disc.visible = visible_now


## Points spread through the water, the same places from frame to frame so that only their number changes.
func _scatter() -> void:
	var random := RandomNumberGenerator.new()
	random.seed = 7
	while _seeds.size() < 2 * MOST_POINTS:
		var p := Vector3(random.randf_range(-1, 1), random.randf_range(-1, 1), random.randf_range(-1, 1))
		if p.length() <= 1.0:
			_seeds.append(p)
	var bottom := -RADIUS + 0.04
	var top := _water_top - 0.03
	for which in 2:
		var many := (_aggregates if which == 0 else _micelles).multimesh
		many.visible_instance_count = _counts[which]
		for i in _counts[which]:
			var s := _seeds[which * MOST_POINTS + i]
			var y := lerpf(bottom, top, 0.5 * (s.y + 1.0))
			var room := sqrt(maxf(RADIUS * RADIUS - y * y, 0.0)) - 0.05
			many.set_instance_transform(i, Transform3D(Basis(), Vector3(s.x * room, y, s.z * room)))


func _show_numbers(frame: RecordingFrame) -> void:
	for child: Node in _numbers.get_children():
		child.queue_free()
	_numbers.add_child(Look.label("At time " + Look.duration(frame.time), 13, Look.MUTED, "bold"))
	var lines := [
		["Oil (and interface)", Look.significant(recording.value(frame, "oil_ml") + recording.value(frame, "interface_ml")) + " mL", Look.AMBER_TEXT],
		["Water phase", Look.significant(recording.value(frame, "water_ml")) + " mL", Look.PRIMARY],
		["Mixed aggregates", Look.moles(recording.value(frame, "mixed_mol")), Look.LIVING_TEXT],
		["Empty micelles", Look.moles(recording.value(frame, "micelles_mol")), Look.VIOLET_TEXT],
		["Hydroxide", Look.moles(recording.value(frame, "hydroxide_mol")), Look.INK],
		["Ethanol", Look.moles(recording.value(frame, "ethanol_mol")), Look.INK]]
	for line: Array in lines:
		var row := HBoxContainer.new()
		var name := Look.label(str(line[0]), 15, Look.MUTED, "medium")
		name.size_flags_horizontal = Control.SIZE_EXPAND_FILL
		row.add_child(name)
		row.add_child(Look.label(str(line[1]), 15, line[2] as Color, "bold"))
		_numbers.add_child(row)
	var fed := recording.value(frame, "inflow_ester")
	var verdict := ""
	if fed > 0.0:
		verdict = "Continuous reactor, residence %s min: %s %% of the ester fed hydrolysed" % [
			Look.significant(recording.value(frame, "residence_min")), Look.significant(100.0 * recording.value(frame, "conversion"))]
	else:
		var cleared := recording.value(frame, "cleared_min")
		verdict = "Clear at %s min (Buhse 1998: 95 min)" % Look.significant(cleared, 4) if not is_nan(cleared) else "Not clear yet"
	var last := Look.label(verdict, 15, Look.LIVING_TEXT, "bold")
	last.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	_numbers.add_child(last)
	_numbers.add_child(Look.label("Balances: %s mol of drift · energy: not closed" % Look.scientific(recording.value(frame, "drift_mol")),
		12, Look.MUTED))


func _process(delta: float) -> void:
	_bar.rotate_y(delta * 6.0)  # a visual turn, not the 800 rpm

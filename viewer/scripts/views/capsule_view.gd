class_name CapsuleView
extends Node3D
## One protocell of Taneja and Higgs, read from a recording (kind "capsule"): its membrane drawn with the shape the model
## computes (a sphere, or a capsule of radius R and cylinder length L), greener as its catalyst is higher; inside, the
## catalyst (amber), the food (blue) and the waste (grey) as points in proportion to their concentrations; at each
## division the daughter let go drifts away and fades, and a burst flashes the membrane red. Over it: the numbers and
## the curves.

const SCALE := 1.1        # scene units for R0
const POINTS := 900       # points inside at most
const RELEASE_FRAMES := 12

var recording: Recording
var _membrane: MeshInstance3D
var _capsule: CapsuleMesh
var _material: StandardMaterial3D
var _released: MeshInstance3D
var _released_material: StandardMaterial3D
var _inside: Array[MultiMeshInstance3D] = []
var _seeds: Array[Vector3] = []
var _card: VBoxContainer
var _chart: CapsuleChart
var _catalyst0: float = NAN


func _init(source: Recording) -> void:
	recording = source


func _ready() -> void:
	var environment := Environment.new()
	environment.background_mode = Environment.BG_COLOR
	environment.background_color = Look.NIGHT
	environment.ambient_light_source = Environment.AMBIENT_SOURCE_COLOR
	environment.ambient_light_color = Color(0.6, 0.66, 0.8)
	environment.ambient_light_energy = 0.7
	environment.tonemap_mode = Environment.TONE_MAPPER_FILMIC
	environment.glow_enabled = true
	environment.glow_intensity = 0.5
	var world := WorldEnvironment.new()
	world.environment = environment
	add_child(world)
	var sun := DirectionalLight3D.new()
	sun.rotation_degrees = Vector3(-45, 30, 0)
	add_child(sun)
	_capsule = CapsuleMesh.new()
	_capsule.radial_segments = 64
	_capsule.rings = 16
	_membrane = MeshInstance3D.new()
	_membrane.mesh = _capsule
	_membrane.rotation_degrees = Vector3(0, 0, 90)  # lying along x
	_material = _shell_material(Color(Look.LIVING, 0.3))
	_membrane.material_override = _material
	add_child(_membrane)
	_released = MeshInstance3D.new()
	var sphere := SphereMesh.new()
	sphere.radial_segments = 48
	sphere.rings = 24
	_released.mesh = sphere
	_released_material = _shell_material(Color(0.75, 0.8, 0.88, 0.25))
	_released.material_override = _released_material
	_released.visible = false
	add_child(_released)
	var random := RandomNumberGenerator.new()
	random.seed = 7
	while _seeds.size() < POINTS:
		_seeds.append(Vector3(random.randf_range(-1, 1), random.randf_range(-1, 1), random.randf_range(-1, 1)))
	for colour: Color in [Look.AMBER, Color(0.35, 0.6, 1.0), Color(0.62, 0.66, 0.72)]:
		_inside.append(_points(colour, 0.03, POINTS))
	var camera := OrbitCamera.new()
	camera.distance = 9.0
	camera.pitch = 0.4
	camera.near_limit = 3.0
	add_child(camera)
	_build_overlay()


func _shell_material(colour: Color) -> StandardMaterial3D:
	var material := StandardMaterial3D.new()
	material.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA
	material.albedo_color = colour
	material.rim_enabled = true
	material.rim = 0.8
	material.roughness = 0.2
	material.emission_enabled = true
	return material


func _points(colour: Color, radius: float, count: int) -> MultiMeshInstance3D:
	var dot := SphereMesh.new()
	dot.radius = radius
	dot.height = 2.0 * radius
	dot.radial_segments = 8
	dot.rings = 4
	var material := StandardMaterial3D.new()
	material.albedo_color = colour
	material.emission_enabled = true
	material.emission = colour
	material.emission_energy_multiplier = 0.7
	dot.material = material
	var many := MultiMesh.new()
	many.transform_format = MultiMesh.TRANSFORM_3D
	many.mesh = dot
	many.instance_count = count
	many.visible_instance_count = 0
	var instance := MultiMeshInstance3D.new()
	instance.multimesh = many
	add_child(instance)
	return instance


func _build_overlay() -> void:
	var layer := CanvasLayer.new()
	add_child(layer)
	var title := VBoxContainer.new()
	title.position = Vector2(40, 84)
	layer.add_child(title)
	title.add_child(Look.label(recording.title, 30, Color.WHITE, "title"))
	var subtitle := Look.label(recording.subtitle, 15, Look.NIGHT_MUTED)
	subtitle.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	subtitle.custom_minimum_size = Vector2(1060, 0)
	title.add_child(subtitle)
	var legend := Look.label("The membrane has the shape the model computes: a sphere, or a capsule of two half-spheres and a "
		+ "cylinder, greener as the catalyst rises. Inside: catalyst (amber), food (blue), waste (grey), in proportion to "
		+ "their concentrations. At a division the daughter let go drifts away; a burst flashes red.", 12, Look.NIGHT_MUTED)
	legend.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	legend.custom_minimum_size = Vector2(1060, 0)
	title.add_child(legend)
	var card := PanelContainer.new()
	card.add_theme_stylebox_override("panel", Look.panel(Color(1, 1, 1, 0.93), 16, 16))
	card.position = Vector2(40, 600)
	card.custom_minimum_size = Vector2(400, 0)
	layer.add_child(card)
	_card = VBoxContainer.new()
	_card.add_theme_constant_override("separation", 4)
	card.add_child(_card)
	_chart = CapsuleChart.new()
	_chart.position = Vector2(1180, 84)
	layer.add_child(_chart)


func show_frame(index: int) -> void:
	var frame := recording.frame(index)
	if frame == null:
		return
	var radius := recording.value(frame, "radius")
	var length := recording.value(frame, "length")
	var catalyst := recording.value(frame, "catalyst")
	if is_nan(_catalyst0):
		_catalyst0 = maxf(recording.value(recording.frame(0), "catalyst"), 1e-6)
	_capsule.radius = radius * SCALE
	_capsule.height = (length + 2.0 * radius) * SCALE
	var life := clampf(catalyst / _catalyst0, 0.0, 1.0)
	_material.albedo_color = Color(Look.LIVING.lerp(Color(0.75, 0.8, 0.88), 1.0 - life), 0.3)
	_material.emission = Color(1.0, 0.3, 0.2) * (1.0 if recording.value(frame, "event") == 2.0 else 0.0)
	_place_inside(frame, radius, length)
	_show_release(index)
	_show_card(frame)
	_chart.show_frame(recording, frame)


## The points inside, scattered through the capsule (a cylinder and its two caps) in proportion to the concentrations.
func _place_inside(frame: RecordingFrame, radius: float, length: float) -> void:
	var keys: Array[String] = ["catalyst", "food", "waste"]
	var offset := 0
	for k in keys.size():
		var many := _inside[k].multimesh
		var count := clampi(int(recording.value(frame, keys[k]) * 250.0), 0, POINTS / 3)
		for i in count:
			var seed := _seeds[(offset + i) % POINTS]
			var along := seed.x * (0.5 * length + radius)
			var cross := Vector2(seed.y, seed.z)
			if cross.length() > 1.0:
				cross = cross.normalized()
			var reach := radius
			var beyond := absf(along) - 0.5 * length
			if beyond > 0.0:
				reach = sqrt(maxf(radius * radius - beyond * beyond, 0.0))
			var point := Vector3(along, cross.x * reach * 0.9, cross.y * reach * 0.9) * SCALE
			many.set_instance_transform(i, Transform3D(Basis(), point))
		many.visible_instance_count = count
		offset += POINTS / 3


## After a division, the daughter let go: a sphere of the rest of the volume, drifting away and fading.
func _show_release(index: int) -> void:
	_released.visible = false
	for back: int in RELEASE_FRAMES:
		var at := index - back
		if at < 0:
			break
		var frame := recording.frame(at)
		if frame == null or recording.value(frame, "event") != 1.0:
			continue
		var share := recording.value(frame, "share")
		var kept := recording.value(frame, "volume")
		var gone := kept / share * (1.0 - share) if share > 0.0 else 0.0
		var r := pow(gone, 1.0 / 3.0) * SCALE
		var progress := float(back) / RELEASE_FRAMES
		_released.visible = true
		_released.scale = Vector3.ONE * r * 2.0
		_released.position = Vector3(0, -1.2 - 4.0 * progress, 1.5 + 3.0 * progress)
		_released_material.albedo_color = Color(0.75, 0.8, 0.88, 0.3 * (1.0 - progress))
		return


func _show_card(frame: RecordingFrame) -> void:
	for child: Node in _card.get_children():
		child.queue_free()
	_card.add_child(Look.label("The protocell", 17, Look.LIVING_TEXT, "bold"))
	var shape := "a sphere" if recording.value(frame, "length") <= 0.0 else "a capsule"
	var lines := [
		["Volume (V0)", Look.significant(recording.value(frame, "volume"), 4)],
		["Lipid area (A0)", Look.significant(recording.value(frame, "lipid_area"), 4)],
		["Reduced volume", Look.significant(recording.value(frame, "reduced"), 3)],
		["Shape", shape],
		["Catalyst inside", Look.significant(recording.value(frame, "catalyst"), 3) + " M"],
		["Food inside", Look.significant(recording.value(frame, "food"), 3) + " M"],
		["Waste inside", Look.significant(recording.value(frame, "waste"), 3) + " M"],
		["Divisions", str(int(recording.value(frame, "divisions")))],
		["Bursts", str(int(recording.value(frame, "bursts")))]]
	for line: Array in lines:
		var row := HBoxContainer.new()
		var name := Look.label(str(line[0]), 14, Look.MUTED, "medium")
		name.size_flags_horizontal = Control.SIZE_EXPAND_FILL
		row.add_child(name)
		row.add_child(Look.label(str(line[1]), 14, Look.INK, "bold"))
		_card.add_child(row)

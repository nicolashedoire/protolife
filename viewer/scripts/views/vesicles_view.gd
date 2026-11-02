class_name VesiclesView
extends Node3D
## Fatty acid vesicles in three dimensions, read from a recording (kind "vesicles"): each group shown by a few of its
## vesicles, shaped by their recorded area and volume (a prolate spheroid of that reduced volume, a sphere when the
## membrane is stretched round its contents), the membrane redder as its tension nears rupture; the kept cargo as
## orange points (always as many, denser when the vesicle shrinks), ribose as blue points inside and out at their
## concentrations, and the endless exchange of the membrane's fatty acids as small points going and coming. Over it:
## a card per group, the scenario's events and the curves.

const PER_GROUP := 3
const BASE := 0.62          # scene units for a 50 nm radius
const CARGO_POINTS := 60
const EXCHANGE_POINTS := 24

var recording: Recording
var _groups: int = 0
var _shells: Array[MeshInstance3D] = []
var _materials: Array[StandardMaterial3D] = []
var _cargo: MultiMeshInstance3D
var _ribose: MultiMeshInstance3D
var _exchange: MultiMeshInstance3D
var _axes: Array[Vector3] = []     # semi-axes of each shown vesicle
var _centres: Array[Vector3] = []
var _cards: Array[VBoxContainer] = []
var _chart: VesiclesChart
var _events: Label
var _time: float = 0.0
var _seeds: Array[Vector3] = []


func _init(source: Recording) -> void:
	recording = source
	_groups = int(recording.constant("records", 1.0))


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
	var random := RandomNumberGenerator.new()
	random.seed = 5
	while _seeds.size() < 2000:
		var p := Vector3(random.randf_range(-1, 1), random.randf_range(-1, 1), random.randf_range(-1, 1))
		if p.length() <= 1.0:
			_seeds.append(p)
	for g in _groups:
		for i in PER_GROUP:
			var shell := MeshInstance3D.new()
			var sphere := SphereMesh.new()
			sphere.radial_segments = 48
			sphere.rings = 24
			shell.mesh = sphere
			var membrane := StandardMaterial3D.new()
			membrane.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA
			membrane.albedo_color = Color(Look.LIVING, 0.28) if g == 0 else Color(0.75, 0.8, 0.88, 0.22)
			membrane.rim_enabled = true
			membrane.rim = 0.8
			membrane.roughness = 0.2
			membrane.emission_enabled = true
			shell.material_override = membrane
			add_child(shell)
			_shells.append(shell)
			_materials.append(membrane)
			_centres.append(Vector3((g - (_groups - 1) / 2.0) * 4.4, (i - 1) * 0.35, (i - 1) * 2.1))
			_axes.append(Vector3.ONE * BASE)
	_cargo = _points(Look.AMBER, 0.035, CARGO_POINTS * PER_GROUP)
	_ribose = _points(Color(0.35, 0.6, 1.0), 0.022, 1200)
	_exchange = _points(Color(0.85, 0.88, 0.95), 0.018, EXCHANGE_POINTS * PER_GROUP * _groups)
	var camera := OrbitCamera.new()
	camera.distance = 8.5
	camera.pitch = 0.45
	camera.near_limit = 3.0
	add_child(camera)
	_build_overlay()


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
	title.add_child(Look.label(recording.subtitle, 15, Look.NIGHT_MUTED))
	_events = Look.label("", 14, Look.NIGHT_MUTED, "medium")
	title.add_child(_events)
	var legend := Look.label("Shapes: an ellipsoid of the same area and volume as the recorded vesicle, stretched or flattened "
		+ "by its reduced volume (an illustration). Orange dots: the cargo; blue: ribose; white: the fatty acids that "
		+ "leave and come back.", 12, Look.NIGHT_MUTED)
	legend.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	legend.custom_minimum_size = Vector2(1060, 0)
	title.add_child(legend)
	for g in _groups:
		var card := PanelContainer.new()
		card.add_theme_stylebox_override("panel", Look.panel(Color(1, 1, 1, 0.93), 16, 16))
		card.position = Vector2(40 + g * 400, 690)
		card.custom_minimum_size = Vector2(380, 0)
		layer.add_child(card)
		var column := VBoxContainer.new()
		column.add_theme_constant_override("separation", 4)
		card.add_child(column)
		_cards.append(column)
	_chart = VesiclesChart.new()
	_chart.position = Vector2(1180, 84)
	layer.add_child(_chart)


## The semi-axes of an ellipsoid of revolution with the area of a sphere of radius `radius` and the reduced volume
## `reduced` (its volume over that sphere's): elongated (prolate) from 0.65 up, flattened (oblate) below, as the
## shapes of least bending energy go (Seifert, Berndl and Lipowsky 1991) — an illustration of the recorded area and
## volume, not a computed shape. The ratio of the axes is found by bisection.
static func spheroid(radius: float, reduced: float) -> Vector3:
	if reduced >= 0.9999:
		return Vector3.ONE * radius
	var prolate := reduced >= 0.65
	var low := 0.0
	var high := 1.0
	var unit_radius := 1.0
	for i in 60:
		var s := 0.5 * (low + high)  # short axis over long axis
		var e := sqrt(maxf(1.0 - s * s, 1e-12))
		var area: float  # over 4 pi times the long (prolate) or equatorial (oblate) axis squared
		var volume: float  # over 4/3 pi times the same axis cubed
		if prolate:
			area = 0.5 * s * s * (1.0 + asin(e) / (s * e))
			volume = s * s
		else:
			area = 0.5 * (1.0 + (1.0 - e * e) / e * atanh(minf(e, 0.999999)))
			volume = s
		unit_radius = sqrt(area)
		if volume / pow(unit_radius, 3.0) > reduced:
			high = s
		else:
			low = s
	var ratio := 0.5 * (low + high)
	var long_axis := radius / unit_radius
	return Vector3(long_axis, long_axis * ratio, long_axis * ratio) if prolate else Vector3(long_axis, long_axis * ratio, long_axis)


func show_frame(index: int) -> void:
	var frame := recording.frame(index)
	if frame == null:
		return
	_time = frame.time
	var rupture := recording.constant("rupture", 10.0)
	for g in _groups:
		var prefix := "r%d_" % g
		var radius := recording.value(frame, prefix + "radius") / 50.0 * BASE
		var reduced := recording.value(frame, prefix + "reduced")
		var tension := recording.value(frame, prefix + "tension")
		var heat := clampf(tension / rupture, 0.0, 1.0)
		for i in PER_GROUP:
			var at := g * PER_GROUP + i
			_axes[at] = spheroid(radius, reduced)
			_shells[at].scale = _axes[at] * 2.0
			_shells[at].position = _centres[at]
			_shells[at].rotation = Vector3(0.0, 0.7 * i, 0.5 + 0.4 * i)
			_materials[at].emission = Color(1.0, 0.35, 0.2) * heat
	_place_cargo(frame)
	_place_ribose(frame)
	_show_cards(frame)
	_chart.show_frame(recording, frame)
	_events.text = _event_text(frame.time)


## Cargo: the same number of points in each vesicle of a group holding cargo, packed in its present volume.
func _place_cargo(frame: RecordingFrame) -> void:
	var many := _cargo.multimesh
	var count := 0
	for g in _groups:
		if recording.value(frame, "r%d_in_cargo" % g) <= 0.0 or is_nan(recording.value(frame, "r%d_in_cargo" % g)):
			continue
		for i in PER_GROUP:
			var at := g * PER_GROUP + i
			for k in CARGO_POINTS:
				var s := _seeds[(k * 7 + i * 131) % _seeds.size()] * 0.85
				var local := Vector3(s.x * _axes[at].x, s.y * _axes[at].y, s.z * _axes[at].z)
				var placed := _shells[at].transform.basis.orthonormalized() * local + _centres[at]
				many.set_instance_transform(count, Transform3D(Basis(), placed))
				count += 1
	many.visible_instance_count = count


## Ribose: points inside each vesicle and in the medium around, as many as their concentration says (10 mM, about four).
func _place_ribose(frame: RecordingFrame) -> void:
	var many := _ribose.multimesh
	var count := 0
	var outside := recording.value(frame, "out_ribose")
	if is_nan(outside):
		many.visible_instance_count = 0
		return
	for k in int(clampf(outside * 4.0, 0.0, 700.0)):
		var s := _seeds[(k * 13) % _seeds.size()]
		many.set_instance_transform(count, Transform3D(Basis(), Vector3(s.x * 6.0, s.y * 2.6, s.z * 4.0)))
		count += 1
	for g in _groups:
		var inside := recording.value(frame, "r%d_in_ribose" % g)
		for i in PER_GROUP:
			var at := g * PER_GROUP + i
			var volume := _axes[at].x * _axes[at].y * _axes[at].z / pow(BASE, 3.0)
			for k in int(clampf(inside * 0.12 * volume, 0.0, 60.0)):
				if count >= many.instance_count:
					break
				var s := _seeds[(k * 17 + at * 59) % _seeds.size()] * 0.8
				var local := Vector3(s.x * _axes[at].x, s.y * _axes[at].y, s.z * _axes[at].z)
				many.set_instance_transform(count, Transform3D(Basis(), _shells[at].transform.basis.orthonormalized() * local + _centres[at]))
				count += 1
	many.visible_instance_count = count


func _process(_delta: float) -> void:
	# fatty acids leaving and coming back: points going out from each membrane and in again, endlessly
	var many := _exchange.multimesh
	var count := 0
	var now := Time.get_ticks_msec() / 1000.0
	for at in _shells.size():
		for k in EXCHANGE_POINTS:
			var direction := _seeds[(k * 23 + at * 41) % _seeds.size()].normalized()
			var phase := fmod(now * 0.6 + k * 0.137 + at * 0.29, 1.0)
			var out := 1.0 + 0.45 * sin(phase * PI)
			var local := Vector3(direction.x * _axes[at].x, direction.y * _axes[at].y, direction.z * _axes[at].z) * out
			many.set_instance_transform(count, Transform3D(Basis(), _shells[at].transform.basis.orthonormalized() * local + _centres[at]))
			count += 1
	many.visible_instance_count = count


func _show_cards(frame: RecordingFrame) -> void:
	for g in _cards.size():
		var column := _cards[g]
		for child: Node in column.get_children():
			child.queue_free()
		var prefix := "r%d_" % g
		var label := str(recording.values[recording.value_index(prefix + "area")]["label"]).split(" : ")[0]
		column.add_child(Look.label(label, 17, Look.LIVING_TEXT if g == 0 else Look.INK, "bold"))
		var lines := [
			["Volume (start = 1)", Look.significant(recording.value(frame, prefix + "volume"), 4)],
			["Area (start = 1)", Look.significant(recording.value(frame, prefix + "area"), 4)],
			["Tension", Look.significant(recording.value(frame, prefix + "tension"), 3) + " mN/m"],
			["Membrane renewed", Look.significant(recording.value(frame, prefix + "turnovers"), 3) + " times"]]
		for key: String in ["cargo", "ribose", "sucrose"]:
			var v := recording.value(frame, prefix + "in_" + key)
			if not is_nan(v):
				lines.append([{"cargo": "Cargo inside", "ribose": "Ribose inside", "sucrose": "Sucrose inside"}[key],
					Look.significant(v, 3) + " mM"])
		for line: Array in lines:
			var row := HBoxContainer.new()
			var name := Look.label(str(line[0]), 14, Look.MUTED, "medium")
			name.size_flags_horizontal = Control.SIZE_EXPAND_FILL
			row.add_child(name)
			row.add_child(Look.label(str(line[1]), 14, Look.INK, "bold"))
			column.add_child(row)


func _event_text(time: float) -> String:
	var parts: PackedStringArray = []
	for key: String in ["event_1", "event_2"]:
		if recording.constants.has(key):
			var at := recording.constant(key)
			parts.append("%s %s: %s" % ["✓" if time >= at else "at", Look.duration(at), recording.constant_labels.get(key, "")])
	return "   ".join(parts)

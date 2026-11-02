class_name WorldView
extends Node3D
## A world of Origin in three dimensions, read from a recording: one of its fields over the grid as a relief, height
## and colour from the frame's own range (logarithmic for amounts, linear for temperature, water and volume
## fractions), the polyester droplets as golden points over their cells, the light warmer when the water is hot.
## Buttons choose the field; the numbers of the frame sit in a panel.

const SIDE := 8.0     # the grid's side in the scene
const HEIGHT := 1.6   # the relief's full height
const MOST_DROPLETS := 3000
const LOGARITHMIC := ["sugar", "food", "activity", "insoluble", "droplets"]

var recording: Recording
var field: int = 0
var _terrain: MeshInstance3D
var _material: StandardMaterial3D
var _droplets: MultiMeshInstance3D
var _sun: DirectionalLight3D
var _environment: Environment
var _buttons: Array[Button] = []
var _legend: Label
var _numbers: VBoxContainer
var _frame: RecordingFrame
var _whole_run: bool = false  # the relief over the run's whole range, not the frame's
var _ranges: Dictionary = {}  # field index -> Vector2(low, high) over the run, in the relief's own scale
var _scale_buttons: Array[Button] = []


func _init(source: Recording, chosen: String) -> void:
	recording = source
	field = maxi(recording.field_index(chosen), 0)


func _ready() -> void:
	_environment = Environment.new()
	_environment.background_mode = Environment.BG_COLOR
	_environment.background_color = Look.NIGHT
	_environment.ambient_light_source = Environment.AMBIENT_SOURCE_COLOR
	_environment.ambient_light_color = Color(0.6, 0.66, 0.78)
	_environment.ambient_light_energy = 0.55
	_environment.tonemap_mode = Environment.TONE_MAPPER_FILMIC
	var world := WorldEnvironment.new()
	world.environment = _environment
	add_child(world)
	_sun = DirectionalLight3D.new()
	_sun.rotation_degrees = Vector3(-55, 40, 0)
	_sun.shadow_enabled = true
	add_child(_sun)
	var slab := MeshInstance3D.new()
	var box := BoxMesh.new()
	box.size = Vector3(SIDE + 0.2, 0.3, SIDE + 0.2)
	slab.mesh = box
	slab.position = Vector3(0, -0.16, 0)
	var stone := StandardMaterial3D.new()
	stone.albedo_color = Color(0.17, 0.2, 0.26)
	slab.material_override = stone
	add_child(slab)
	_terrain = MeshInstance3D.new()
	_material = StandardMaterial3D.new()
	_material.vertex_color_use_as_albedo = true
	_material.roughness = 0.65
	_terrain.material_override = _material
	add_child(_terrain)
	_droplets = _golden_points()
	var camera := OrbitCamera.new()
	camera.distance = 12.0
	camera.pitch = 0.75
	camera.far_limit = 40.0
	camera.near_limit = 3.0
	add_child(camera)
	_build_overlay()


func _golden_points() -> MultiMeshInstance3D:
	var dot := SphereMesh.new()
	dot.radius = 0.045
	dot.height = 0.09
	dot.radial_segments = 8
	dot.rings = 5
	var gold := StandardMaterial3D.new()
	gold.albedo_color = Look.AMBER
	gold.emission_enabled = true
	gold.emission = Look.AMBER
	gold.emission_energy_multiplier = 0.8
	dot.material = gold
	var many := MultiMesh.new()
	many.transform_format = MultiMesh.TRANSFORM_3D
	many.mesh = dot
	many.instance_count = MOST_DROPLETS
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
	var choices := HBoxContainer.new()
	choices.position = Vector2(40, 170)
	choices.add_theme_constant_override("separation", 8)
	layer.add_child(choices)
	for i in recording.fields.size():
		var button := Look.button(str(recording.fields[i]["label"]), i == field)
		button.pressed.connect(_choose.bind(i))
		_buttons.append(button)
		choices.add_child(button)
	var scales := HBoxContainer.new()
	scales.position = Vector2(40, 214)
	scales.add_theme_constant_override("separation", 8)
	layer.add_child(scales)
	for i in 2:
		var button := Look.button(["Relief of this frame", "Scale of the whole run"][i], (i == 1) == _whole_run)
		button.pressed.connect(_scale.bind(i == 1))
		_scale_buttons.append(button)
		scales.add_child(button)
	var panel := PanelContainer.new()
	panel.add_theme_stylebox_override("panel", Look.panel(Color(1, 1, 1, 0.93), 16, 18))
	panel.position = Vector2(1440, 84)
	panel.custom_minimum_size = Vector2(440, 0)
	layer.add_child(panel)
	var column := VBoxContainer.new()
	column.add_theme_constant_override("separation", 6)
	panel.add_child(column)
	_legend = Look.label("", 14, Look.INK, "medium")
	_legend.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	_legend.custom_minimum_size = Vector2(400, 0)
	column.add_child(_legend)
	_numbers = VBoxContainer.new()
	column.add_child(_numbers)


func _choose(which: int) -> void:
	field = which
	for i in _buttons.size():
		_restyle(_buttons[i], i == field)
	if _frame != null:
		_draw_field()


func _scale(whole: bool) -> void:
	_whole_run = whole
	for i in _scale_buttons.size():
		_restyle(_scale_buttons[i], (i == 1) == _whole_run)
	if _frame != null:
		_draw_field()


func _restyle(button: Button, active: bool) -> void:
	var fresh := Look.button(button.text, active)
	for style: String in ["normal", "hover", "pressed"]:
		button.add_theme_stylebox_override(style, fresh.get_theme_stylebox(style))
	for colour: String in ["font_color", "font_hover_color", "font_pressed_color"]:
		button.add_theme_color_override(colour, fresh.get_theme_color(colour))
	fresh.free()


## A field in the relief's scale: the decimal logarithm of an amount (its smallest positive value standing for zero),
## or the value itself.
func _scaled(raw: PackedFloat32Array, logarithmic: bool) -> PackedFloat32Array:
	var shown := PackedFloat32Array()
	shown.resize(raw.size())
	var positive := INF
	for v in raw:
		if v > 0.0:
			positive = minf(positive, v)
	for c in raw.size():
		shown[c] = log(maxf(raw[c], positive if positive < INF else 1e-30)) / log(10.0) if logarithmic else raw[c]
	return shown


## The field's range over every frame of the run, in the relief's scale (read once).
func _run_range(which: int, logarithmic: bool) -> Vector2:
	if _ranges.has(which):
		return _ranges[which]
	var low := INF
	var high := -INF
	var cells := recording.n * recording.n
	for i in recording.frame_count():
		for v in _scaled(recording.frame(i).field(which, cells), logarithmic):
			low = minf(low, v)
			high = maxf(high, v)
	_ranges[which] = Vector2(low, high)
	return _ranges[which]


func show_frame(index: int) -> void:
	_frame = recording.frame(index)
	if _frame == null:
		return
	_draw_field()
	_draw_droplets()
	_show_numbers()
	var kelvin := recording.value(_frame, "mean_temperature")
	var warmth := clampf((kelvin - 293.15) / 60.0, 0.0, 1.0) if not is_nan(kelvin) else 0.0
	_sun.light_color = Color(0.85, 0.9, 1.0).lerp(Color(1.0, 0.78, 0.5), warmth)
	_sun.light_energy = 0.9 + 0.6 * warmth


## Heights and colours of the chosen field over the grid, from the frame's own range.
func _draw_field() -> void:
	var n := recording.n
	var cells := n * n
	var key := str(recording.fields[field]["key"])
	var logarithmic := key in LOGARITHMIC
	var shown := _scaled(_frame.field(field, cells), logarithmic)
	var low := INF
	var high := -INF
	for v in shown:
		low = minf(low, v)
		high = maxf(high, v)
	if _whole_run:
		var whole := _run_range(field, logarithmic)
		low = whole.x
		high = whole.y
	var span := high - low
	var tool := SurfaceTool.new()
	tool.begin(Mesh.PRIMITIVE_TRIANGLES)
	var step := SIDE / maxf(n - 1, 1)
	for row in n:
		for col in n:
			var t := (shown[row * n + col] - low) / span if span > 0.0 else 0.0
			tool.set_color(_colour(t))
			tool.add_vertex(Vector3(-SIDE / 2 + col * step, 0.02 + HEIGHT * t, -SIDE / 2 + row * step))
	for row in n - 1:
		for col in n - 1:
			var a := row * n + col
			tool.add_index(a)
			tool.add_index(a + 1)
			tool.add_index(a + n)
			tool.add_index(a + 1)
			tool.add_index(a + n + 1)
			tool.add_index(a + n)
	tool.generate_normals()
	_terrain.mesh = tool.commit()
	var unit := str(recording.fields[field]["unit"])
	var lowest := _undo(low, logarithmic)
	var highest := _undo(high, logarithmic)
	_legend.text = "%s (%s), from %s to %s %s; %s, %s" % [recording.fields[field]["label"], unit, Look.significant(lowest, 4),
		Look.significant(highest, 4), "over the whole run" if _whole_run else "in this frame",
		"logarithmic scale" if logarithmic else "linear scale",
		"a flat relief: a uniform field" if high - low <= 0.0 else "the height covers this range"]


func _undo(v: float, logarithmic: bool) -> float:
	return pow(10.0, v) if logarithmic and v > -INF and v < INF else v


## Deep blue to teal to warm yellow.
func _colour(t: float) -> Color:
	var stops := [Color(0.12, 0.18, 0.45), Color(0.06, 0.62, 0.62), Color(0.55, 0.82, 0.35), Color(0.99, 0.86, 0.25)]
	var scaled := clampf(t, 0.0, 1.0) * (stops.size() - 1)
	var i := mini(int(scaled), stops.size() - 2)
	return (stops[i] as Color).lerp(stops[i + 1] as Color, scaled - i)


## The droplets per mL of each cell, as points over it (a few per cell, more where they are dense).
func _draw_droplets() -> void:
	var which := recording.field_index("droplets")
	var many := _droplets.multimesh
	if which < 0:
		many.visible_instance_count = 0
		return
	var n := recording.n
	var density := _frame.field(which, n * n)
	var densest := 0.0
	for d in density:
		densest = maxf(densest, d)
	var random := RandomNumberGenerator.new()
	random.seed = 11
	var count := 0
	var step := SIDE / maxf(n - 1, 1)
	for c in n * n:
		if densest <= 0.0 or density[c] <= 0.0:
			continue
		var here := 1 + int(3.0 * density[c] / densest)
		for k in here:
			if count >= MOST_DROPLETS:
				break
			var spot := Vector3(-SIDE / 2 + (c % n) * step + random.randf_range(-0.4, 0.4) * step,
				HEIGHT + 0.25 + random.randf_range(0.0, 0.3), -SIDE / 2 + (c / n) * step + random.randf_range(-0.4, 0.4) * step)
			many.set_instance_transform(count, Transform3D(Basis(), spot))
			count += 1
	many.visible_instance_count = count


func _show_numbers() -> void:
	for child: Node in _numbers.get_children():
		child.queue_free()
	_numbers.add_child(Look.label("At time " + Look.duration(_frame.time), 13, Look.MUTED, "bold"))
	for i in recording.values.size():
		var quantity := recording.values[i]
		var row := HBoxContainer.new()
		var name := Look.label(str(quantity["label"]), 14, Look.MUTED, "medium")
		name.size_flags_horizontal = Control.SIZE_EXPAND_FILL
		row.add_child(name)
		row.add_child(Look.label(Look.significant(_frame.values[i]) + " " + str(quantity["unit"]), 14, Look.INK, "bold"))
		_numbers.add_child(row)

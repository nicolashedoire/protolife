class_name LineageView
extends Node3D
## Vesicles that grow on micelles and divide under shear, read from a recording (kind "lineage", entries 75 and 76; and
## entries 83 to 85, where the inside is made by Origin's chemistry: the constant inside_sugars, the cargo being sugars).
## FOUNDERS vesicles of the start and all their descendants are shown, each daughter on a ring round the place of its
## mother, so that the family shows as nested rings. Each vesicle is shaped by its recorded area and reduced volume.
## The founders' cargo is shown as orange points shared out among the daughters: a vesicle holding a fifth of a
## founder's cargo shows a fifth of its points, the last one smaller in proportion. Vesicles born of the micelles,
## empty, are grey and around (with a chemistry, tinted amber as far as their inside reaches the lineage's); micelles
## waiting in the medium are pale points; the shear is drawn as streaks sliding faster the higher they are. Over it: a
## card each for the lineage, the vesicles born of the micelles and the medium; the feeds and shears; the curves.

const FOUNDERS := 4
const BASE := 0.62          # scene units for a 50 nm radius
const CARGO_POINTS := 40    # per vesicle of the start
const MOST := 1600          # vesicles shown at most
const STREAKS := 160

var recording: Recording
var _groups: int = 0
var _sugars: bool = false           # the inside is made by a chemistry (entries 83 to 85), not a kept cargo
var _history: Array[Array] = []     # per group, each division: Vector2(daughters, radius after in nm)
var _scanned: int = 0
var _last_generation := PackedInt32Array()
var _last_count := PackedFloat32Array()
var _widest: float = 50.0           # nm, the largest radius the lineage reaches
var _spacing: float = 4.0
var _membranes: MultiMeshInstance3D
var _cargo: MultiMeshInstance3D
var _micelles: MultiMeshInstance3D
var _streaks: MultiMeshInstance3D
var _streak_x := PackedFloat32Array()
var _seeds: Array[Vector3] = []
var _shear: float = 0.0
var _cards: Array[VBoxContainer] = []
var _chart: LineageChart
var _events: Label


func _init(source: Recording) -> void:
	recording = source
	_groups = int(recording.constant("records", 1.0))
	_sugars = recording.constant("inside_sugars", 0.0) > 0.5
	for g in _groups:
		_history.append([])
		_last_generation.append(-1)
		_last_count.append(0.0)


func _ready() -> void:
	var environment := Environment.new()
	environment.background_mode = Environment.BG_COLOR
	environment.background_color = Look.NIGHT
	environment.ambient_light_source = Environment.AMBIENT_SOURCE_COLOR
	environment.ambient_light_color = Color(0.6, 0.66, 0.8)
	environment.ambient_light_energy = 0.7
	environment.tonemap_mode = Environment.TONE_MAPPER_FILMIC
	environment.glow_enabled = true
	environment.glow_intensity = 0.4
	var world := WorldEnvironment.new()
	world.environment = environment
	add_child(world)
	var sun := DirectionalLight3D.new()
	sun.rotation_degrees = Vector3(-50, 30, 0)
	add_child(sun)
	var random := RandomNumberGenerator.new()
	random.seed = 7
	while _seeds.size() < 3000:
		var p := Vector3(random.randf_range(-1, 1), random.randf_range(-1, 1), random.randf_range(-1, 1))
		if p.length() <= 1.0:
			_seeds.append(p)
	_scan()
	_spacing = _extent() * 2.4
	_membranes = _many(_membrane_mesh(), MOST, true)
	_cargo = _many(_dot(Look.AMBER, 0.05), MOST, false)
	_micelles = _many(_dot(Color(0.95, 0.93, 0.8), 0.03), 900, false)
	var streak := BoxMesh.new()
	streak.size = Vector3(0.9, 0.012, 0.012)
	var faint := StandardMaterial3D.new()
	faint.shading_mode = BaseMaterial3D.SHADING_MODE_UNSHADED
	faint.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA
	faint.albedo_color = Color(Look.VIOLET, 0.55)
	streak.material = faint
	_streaks = _many(streak, STREAKS, false)
	for i in STREAKS:
		_streak_x.append(_seeds[i * 11 % _seeds.size()].x)
	var camera := OrbitCamera.new()
	camera.distance = _spacing * (5.0 if _groups > 1 else 4.0)
	camera.far_limit = _spacing * 8.0
	camera.near_limit = _spacing * 0.6
	camera.pitch = 0.6
	add_child(camera)
	_build_overlay()


func _membrane_mesh() -> SphereMesh:
	var sphere := SphereMesh.new()
	sphere.radial_segments = 32
	sphere.rings = 16
	var membrane := StandardMaterial3D.new()
	membrane.transparency = BaseMaterial3D.TRANSPARENCY_ALPHA
	membrane.vertex_color_use_as_albedo = true
	membrane.rim_enabled = true
	membrane.rim = 0.8
	membrane.roughness = 0.2
	sphere.material = membrane
	return sphere


func _dot(colour: Color, radius: float) -> SphereMesh:
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
	return dot


func _many(mesh: Mesh, count: int, colours: bool) -> MultiMeshInstance3D:
	var many := MultiMesh.new()
	many.transform_format = MultiMesh.TRANSFORM_3D
	many.use_colors = colours
	many.mesh = mesh
	many.instance_count = count
	many.visible_instance_count = 0
	var instance := MultiMeshInstance3D.new()
	instance.multimesh = many
	add_child(instance)
	return instance


## Reads the frames not read yet for each group's divisions: when its generation goes up, how many daughters each
## vesicle made (the jump in its count) and their radius.
func _scan() -> void:
	var count := recording.frame_count()
	while _scanned < count:
		var frame := recording.frame(_scanned)
		_scanned += 1
		for g in _groups:
			var prefix := "r%d_" % g
			if recording.value(frame, prefix + "present") < 0.5:
				continue
			if g == 0:
				_widest = maxf(_widest, recording.value(frame, prefix + "radius"))
			var generation := int(recording.value(frame, prefix + "generation"))
			var now := recording.value(frame, prefix + "count")
			if _last_generation[g] >= 0 and generation > _last_generation[g]:
				var steps := generation - _last_generation[g]
				var daughters := float(roundi(pow(now / _last_count[g], 1.0 / steps)))
				for k in steps:
					_history[g].append(Vector2(daughters, recording.value(frame, prefix + "radius")))
			_last_generation[g] = generation
			_last_count[g] = now


## The ring on which a division puts its daughters, round their mother's place.
static func ring(daughters: float, radius_nm: float) -> float:
	return 0.0 if daughters < 2.0 else radius_nm / 50.0 * BASE * 1.3 / sin(PI / daughters)


## How far the lineage of one vesicle of the start reaches from its place, over all its divisions.
func _extent() -> float:
	var reach := _widest / 50.0 * BASE * 1.6  # a flattened vesicle's long axis, about
	for step: Vector2 in _history[0]:
		reach += ring(step.x, step.y)
	return reach


func _lineage_at() -> Vector3:
	return Vector3(-0.9 * _spacing, 0.0, 0.0) if _groups > 1 else Vector3.ZERO


## The places of the lineage's vesicles at a generation: the founders, then each division's daughters on a ring round
## their mother's place, the ring turned a quarter at each generation so that the family shows in three dimensions.
func _places(generation: int) -> PackedVector3Array:
	var places := PackedVector3Array()
	for f in FOUNDERS:  # a square of founders, beside the cloud of vesicles born empty when there is one
		places.append(_lineage_at() + Vector3((f % 2) - 0.5, 0.0, (f / 2) - 0.5) * _spacing)
	var history: Array = _history[0]
	for level in mini(generation, history.size()):
		var step: Vector2 = history[level]
		var daughters := int(step.x)
		if daughters < 2:
			continue
		var radius := ring(step.x, step.y)
		var next := PackedVector3Array()
		for p: Vector3 in places:
			for j in daughters:
				var angle := TAU * j / daughters + level * 0.7
				var offset := Vector3(cos(angle), 0.0, sin(angle)) if level % 2 == 0 else Vector3(cos(angle), sin(angle), 0.0)
				next.append(p + offset * radius)
		places = next
	return places


func _shape(prefix: String, frame: RecordingFrame) -> Vector3:
	return VesiclesView.spheroid(recording.value(frame, prefix + "radius") / 50.0 * BASE, recording.value(frame, prefix + "reduced"))


func _turn(index: int) -> Basis:
	var s := _seeds[(index * 31 + 3) % _seeds.size()]
	return Basis.from_euler(s * PI)


func show_frame(index: int) -> void:
	_scan()
	var frame := recording.frame(index)
	if frame == null:
		return
	var rupture := recording.constant("rupture", 10.0)
	var membranes := _membranes.multimesh
	var cargo := _cargo.multimesh
	var shown := 0
	var points := 0
	var axes := _shape("r0_", frame)
	var heat := clampf(recording.value(frame, "r0_tension") / rupture, 0.0, 1.0)
	var places := _places(int(recording.value(frame, "r0_generation")))
	var per_vesicle := CARGO_POINTS * recording.value(frame, "r0_cargo_share")
	var dot := maxf(1.0, recording.value(frame, "r0_radius") / 50.0)  # points grow with vesicles past 50 nm
	for i in mini(places.size(), MOST):
		var turn := _turn(i)
		membranes.set_instance_transform(shown, Transform3D(turn * Basis.from_scale(axes * 2.0), places[i]))
		membranes.set_instance_color(shown, Color(Look.LIVING.lerp(Look.DANGER, heat), 0.3))
		shown += 1
		if is_nan(per_vesicle):
			continue
		var whole := int(floor(per_vesicle + 1e-6))
		var rest := per_vesicle - whole
		for k in whole + (1 if rest > 0.02 else 0):
			if points >= MOST:
				break
			var s := _seeds[(k * 7 + i * 131) % _seeds.size()] * 0.8
			var size := dot * (1.0 if k < whole else pow(rest, 1.0 / 3.0))
			cargo.set_instance_transform(points, Transform3D(Basis.from_scale(Vector3.ONE * size),
				places[i] + turn * Vector3(s.x * axes.x, s.y * axes.y, s.z * axes.z)))
			points += 1
	shown = _place_born(frame, shown)
	membranes.visible_instance_count = shown
	cargo.visible_instance_count = points
	_place_micelles(frame)
	_shear = recording.value(frame, "shear")
	_show_cards(frame)
	_chart.show_frame(recording, frame)
	_events.text = _event_text(frame.time)


## The other vesicles (born of the micelles; with a chemistry, also the fresh ones of the start): as many as their groups
## count, in a cloud beside the lineage; grey when empty, amber as far as their inside reaches the lineage's.
func _place_born(frame: RecordingFrame, shown: int) -> int:
	var membranes := _membranes.multimesh
	var at := shown
	var lineage := recording.value(frame, "r0_cargo")
	for g in range(1, _groups):
		var prefix := "r%d_" % g
		if recording.value(frame, prefix + "present") < 0.5:
			continue
		var axes := _shape(prefix, frame)
		for i in roundi(recording.value(frame, prefix + "count") * FOUNDERS):
			if at >= MOST:
				return at
			var s := _seeds[(g * 977 + i * 13) % _seeds.size()]
			var place := Vector3(1.5 * _spacing, 0.0, 0.0) + Vector3(s.x, s.y * 0.5, s.z) * _spacing
			membranes.set_instance_transform(at, Transform3D(_turn(at) * Basis.from_scale(axes * 2.0), place))
			var grey := Color(0.75, 0.8, 0.88, 0.22)
			var lit := clampf(recording.value(frame, prefix + "cargo") / lineage, 0.0, 1.0) if _sugars and lineage > 0.0 else 0.0
			membranes.set_instance_color(at, grey.lerp(Color(Look.AMBER, 0.4), 0.0 if is_nan(lit) else lit))
			at += 1
	return at


## Micelles waiting in the medium (on the membranes, as intermediate, or making vesicles): about fifty points a mM.
func _place_micelles(frame: RecordingFrame) -> void:
	var waiting := recording.value(frame, "coat") + recording.value(frame, "intermediate") + recording.value(frame, "forming")
	var many := _micelles.multimesh
	var count := int(clampf(waiting * 50.0, 0.0, many.instance_count))
	for k in count:
		var s := _seeds[(k * 17 + 5) % _seeds.size()]
		many.set_instance_transform(k, Transform3D(Basis.from_scale(Vector3.ONE * _grain()),
			Vector3(s.x * 3.2, s.y * 1.4, s.z * 3.2) * _spacing))
	many.visible_instance_count = count


## Small marks (micelles, streaks) grow with a scene wider than entry 76's.
func _grain() -> float:
	return maxf(1.0, _spacing / 10.0)


func _process(delta: float) -> void:
	# the shear: streaks sliding along x, faster the higher they are, as the layers of a sheared liquid
	var many := _streaks.multimesh
	if _shear <= 0.0:
		many.visible_instance_count = 0
		return
	for i in STREAKS:
		var s := _seeds[(i * 11) % _seeds.size()]
		var height := s.y * 1.6 * _spacing
		_streak_x[i] = fmod(_streak_x[i] + delta * 0.12 * s.y * _shear / 20.0 + 3.0, 2.0) - 1.0
		many.set_instance_transform(i, Transform3D(Basis.from_scale(Vector3.ONE * _grain()),
			Vector3(_streak_x[i] * 3.4 * _spacing, height, s.z * 3.0 * _spacing)))
	many.visible_instance_count = STREAKS


func _build_overlay() -> void:
	var layer := CanvasLayer.new()
	add_child(layer)
	var title := VBoxContainer.new()
	title.position = Vector2(40, 84)
	layer.add_child(title)
	for line: Label in [Look.label(recording.title, 30, Color.WHITE, "title"), Look.label(recording.subtitle, 15, Look.NIGHT_MUTED)]:
		line.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART  # the column stops short of the chart
		line.custom_minimum_size = Vector2(1060, 0)
		title.add_child(line)
	_events = Look.label("", 14, Look.NIGHT_MUTED, "medium")
	title.add_child(_events)
	var text := ("%d starting vesicles and all their descendants, each daughter on a ring around its mother; "
		+ "shapes: an ellipsoid of the same area and volume (an illustration). Orange dots: the cargo, %d per starting "
		+ "vesicle, shared between the daughters (the last, smaller dot carries the rest). Grey: vesicles born empty of the "
		+ "micelles. Pale dots: micelles waiting. Purple strokes: the shear.") % [FOUNDERS, CARGO_POINTS]
	if _sugars:
		text = ("%d seeded vesicles and all their descendants, each daughter on a ring around its mother; shapes: an "
			+ "ellipsoid of the same area and volume (an illustration); red membrane: stretched to rupture. Orange "
			+ "dots: the sugars the chemistry makes, %d per seeded vesicle at the first cycle; shared at division, made again "
			+ "afterwards if there is a chemistry. On the right, the other vesicles: grey if empty, orange if their chemistry lights up. "
			+ "Pale dots: micelles waiting. Purple strokes: the shear.") % [FOUNDERS, CARGO_POINTS]
	var legend := Look.label(text, 12, Look.NIGHT_MUTED)
	legend.autowrap_mode = TextServer.AUTOWRAP_WORD_SMART
	legend.custom_minimum_size = Vector2(1060, 0)
	title.add_child(legend)
	for c in 3:
		var card := PanelContainer.new()
		card.add_theme_stylebox_override("panel", Look.panel(Color(1, 1, 1, 0.93), 16, 16))
		card.position = Vector2(40 + c * 400, 760)
		card.custom_minimum_size = Vector2(380, 0)
		layer.add_child(card)
		var column := VBoxContainer.new()
		column.add_theme_constant_override("separation", 4)
		card.add_child(column)
		_cards.append(column)
	_chart = LineageChart.new()
	_chart.position = Vector2(1180, 84)
	layer.add_child(_chart)


func _show_cards(frame: RecordingFrame) -> void:
	var share := recording.value(frame, "r0_cargo_share")
	var born_share := 0.0
	var groups := 0
	for g in range(1, _groups):
		if recording.value(frame, "r%d_present" % g) > 0.5 and recording.value(frame, "r%d_origin" % g) > 0.5:
			born_share += recording.value(frame, "r%d_share" % g)
			groups += 1
	var waiting := recording.value(frame, "coat") + recording.value(frame, "intermediate") + recording.value(frame, "forming")
	var held := "—"
	if not is_nan(share):
		if _sugars:
			held = Look.significant(share, 3) + " of a seeded one"
		else:
			held = "1/%d of the start" % roundi(1.0 / share) if share < 0.999 else "the starting one"
	var inside := "—" if is_nan(share) else Look.significant(recording.value(frame, "r0_cargo"), 3) + (" M" if _sugars else " mM")
	var cards := [
		["The lineage of the seeded vesicles" if _sugars else "The lineage of the starting vesicles", Look.LIVING_TEXT, [
			["Descendants per starting vesicle", Look.significant(recording.value(frame, "r0_count"), 4)],
			["Generation", str(int(recording.value(frame, "r0_generation")))],
			["Radius", Look.significant(recording.value(frame, "r0_radius"), 3) + " nm"],
			["Reduced volume", Look.significant(recording.value(frame, "r0_reduced"), 3)],
			["Sugars per vesicle" if _sugars else "Cargo per vesicle", held],
			["Sugar carbon inside" if _sugars else "Cargo inside", inside],
			["Share of the membranes", Look.significant(100.0 * recording.value(frame, "r0_share"), 3) + " %"]]],
		["Born of micelles" if _sugars else "Born of micelles, empty", Look.INK, _born_lines(frame, groups, born_share)],
		["The medium", Look.INK, [
			["Micelles waiting", (Look.significant(waiting, 3) if waiting > 1e-6 else "0") + " mM"],
			["Fatty acid added", Look.significant(recording.value(frame, "fed"), 3) + " equivalent(s)"],
			["Shear", Look.significant(recording.value(frame, "shear"), 3) + " s⁻¹"],
			["Vesicles born of a division", Look.significant(recording.value(frame, "divided"), 3) + " per starting vesicle"]]]]
	if recording.value_index("medium_sugar") >= 0:
		(cards[2][2] as Array).push_front(["Sugar carbon in the medium", Look.significant(recording.value(frame, "medium_sugar"), 3) + " M"])
	for c in _cards.size():
		var column := _cards[c]
		for child: Node in column.get_children():
			child.queue_free()
		var card: Array = cards[c]
		column.add_child(Look.label(str(card[0]), 17, card[1], "bold"))
		for line: Array in card[2]:
			var row := HBoxContainer.new()
			var name := Look.label(str(line[0]), 14, Look.MUTED, "medium")
			name.size_flags_horizontal = Control.SIZE_EXPAND_FILL
			row.add_child(name)
			row.add_child(Look.label(str(line[1]), 14, Look.INK, "bold"))
			column.add_child(row)


func _born_lines(frame: RecordingFrame, groups: int, born_share: float) -> Array:
	var lines := [["Vesicles per starting vesicle", Look.significant(LineageChart.born(recording, frame), 3)],
		["Feeds that gave birth", str(groups)],
		["Share of the membranes", Look.significant(100.0 * born_share, 3) + " %"]]
	if _sugars:
		var carbon := 0.0
		var weight := 0.0
		for g in range(1, _groups):
			var prefix := "r%d_" % g
			if recording.value(frame, prefix + "present") > 0.5 and recording.value(frame, prefix + "origin") > 0.5:
				carbon += recording.value(frame, prefix + "share") * recording.value(frame, prefix + "cargo")
				weight += recording.value(frame, prefix + "share")
		lines[1] = ["Sugar carbon inside", Look.significant(carbon / weight, 3) + " M" if weight > 0.0 else "—"]
	return lines


func _event_text(time: float) -> String:
	var feeds := int(recording.constant("feeds", 1.0))
	var doses := "%d doses %s apart" % [feeds, Look.duration(recording.constant("interval"))] if feeds > 1 else "one dose"
	if not recording.constants.has("cycles"):
		return "%s equivalent(s) added at %s, in %s" % [Look.significant(recording.constant("micelles"), 3),
			Look.duration(recording.constant("first_feed")), doses]
	var period := recording.constant("period")
	var cycles := int(recording.constant("cycles"))
	var cycle := mini(int(time / period) + 1, cycles)
	var shear_at := recording.constant("shear_at")
	if recording.constant("stirred", 0.0) > 0.5:
		return "Cycle %d of %d · food %s into each cycle, in %s · shear of %s all along" % [cycle, cycles,
			Look.duration(recording.constant("first_feed")), doses, Look.significant(recording.constant("shear_rate"), 3) + " s⁻¹"]
	var text := "Cycle %d of %d · food from %s, in %s · shear of %s for %s, %s into each cycle" % [cycle, cycles,
		Look.duration(recording.constant("first_feed")), doses, Look.significant(recording.constant("shear_rate"), 3) + " s⁻¹",
		Look.duration(recording.constant("shear_for")), Look.duration(shear_at)]
	if _shear > 0.0:
		text += "   ▶ shear under way"
	return text

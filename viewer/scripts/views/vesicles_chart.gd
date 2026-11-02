class_name VesiclesChart
extends Control
## The vesicles over time, up to the frame shown: each group's volume and area against the start, and the cargo's
## concentration against its start where there is one; the scenario's events as vertical marks.

const COLOURS := [Color8(16, 179, 154), Color8(120, 132, 150)]

var _times := PackedFloat32Array()
var _series: Array[PackedFloat32Array] = []
var _names: PackedStringArray = []
var _colours: Array[Color] = []
var _dashed: Array[bool] = []
var _keys: PackedStringArray = []
var _read: int = 0
var _current: int = 0
var _span: float = 60.0
var _events: PackedFloat32Array = []
var _cargo0: float = NAN


func _init() -> void:
	custom_minimum_size = Vector2(700, 300)


func _prepare(recording: Recording) -> void:
	var groups := int(recording.constant("records", 1.0))
	for g in groups:
		var group := str(recording.values[recording.value_index("r%d_volume" % g)]["label"]).split(" : ")[0]
		_add("r%d_volume" % g, group + ": volume", COLOURS[mini(g, 1)], false)
		_add("r%d_area" % g, group + ": area", COLOURS[mini(g, 1)], true)
	if recording.value_index("r0_in_cargo") >= 0:
		_add("r0_in_cargo", "Cargo (concentration)", Look.AMBER, false)
	for key: String in ["event_1", "event_2"]:
		if recording.constants.has(key):
			_events.append(recording.constant(key))


func _add(key: String, name: String, colour: Color, dashed: bool) -> void:
	_keys.append(key)
	_names.append(name)
	_colours.append(colour)
	_dashed.append(dashed)
	_series.append(PackedFloat32Array())


func show_frame(recording: Recording, frame: RecordingFrame) -> void:
	if _keys.is_empty():
		_prepare(recording)
	var count := recording.frame_count()
	while _read < count:
		var f := recording.frame(_read)
		_read += 1
		_times.append(f.time)
		for i in _keys.size():
			var v := recording.value(f, _keys[i])
			if _keys[i] == "r0_in_cargo":
				if is_nan(_cargo0):
					_cargo0 = v
				v = v / _cargo0 if _cargo0 > 0.0 else NAN
			_series[i].append(v)
	_current = frame.index
	_span = maxf(60.0, _times[_times.size() - 1] if not _times.is_empty() else 60.0)
	queue_redraw()


func _draw() -> void:
	draw_style_box(Look.panel(Color(1, 1, 1, 0.93), 16, 0), Rect2(Vector2.ZERO, size))
	draw_string(Look.font("title"), Vector2(20, 30), "Over time (start = 1)", HORIZONTAL_ALIGNMENT_LEFT, -1, 17, Look.INK)
	var x := 20.0
	var y := 52.0
	for i in _names.size():
		var width := Look.font("medium").get_string_size(_names[i], HORIZONTAL_ALIGNMENT_LEFT, -1, 11).x + 34.0
		if x + width > size.x - 20:
			x = 20.0
			y += 16.0
		draw_rect(Rect2(x, y - 5, 16, 3 if not _dashed[i] else 2), _colours[i])
		draw_string(Look.font("medium"), Vector2(x + 22, y), _names[i], HORIZONTAL_ALIGNMENT_LEFT, -1, 11, Look.INK)
		x += width
	var plot := Rect2(54, y + 16, size.x - 80, size.y - y - 50)
	draw_rect(plot, Look.LINE, false, 1.0)
	var low := 0.6
	var high := 1.45
	for mark: float in [0.75, 1.0, 1.25]:
		var at_y := plot.end.y - (mark - low) / (high - low) * plot.size.y
		draw_line(Vector2(plot.position.x, at_y), Vector2(plot.end.x, at_y), Look.LINE, 1.0)
		draw_string(Look.font("mono"), Vector2(10, at_y + 4), Look.significant(mark, 3), HORIZONTAL_ALIGNMENT_LEFT, -1, 11, Look.MUTED)
	for tick: Dictionary in Look.ticks(_span):
		var at_x := plot.position.x + float(tick["time"]) / _span * plot.size.x
		draw_string(Look.font("mono"), Vector2(at_x - 20, plot.end.y + 16), str(tick["text"]), HORIZONTAL_ALIGNMENT_CENTER, 40, 11,
			Look.MUTED)
	for event in _events:
		var at := plot.position.x + event / _span * plot.size.x
		draw_line(Vector2(at, plot.position.y), Vector2(at, plot.end.y), Look.AMBER, 1.2)
	for i in _series.size():
		var points := PackedVector2Array()
		for k in mini(_current + 1, _times.size()):
			var v := _series[i][k]
			if is_nan(v):
				continue
			points.append(Vector2(plot.position.x + _times[k] / _span * plot.size.x,
				plot.end.y - (clampf(v, low, high) - low) / (high - low) * plot.size.y))
		if points.size() > 1:
			draw_polyline(points, _colours[i], 1.4 if _dashed[i] else 2.4, true)

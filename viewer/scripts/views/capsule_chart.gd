class_name CapsuleChart
extends Control
## The protocell over time, up to the frame shown: its volume and lipid area (logarithmic scale, against V0 and A0),
## its reduced volume, and its catalyst; each division as a vertical mark.

var _times := PackedFloat32Array()
var _series: Array[PackedFloat32Array] = []
var _divisions := PackedFloat32Array()
var _read: int = 0
var _current: int = 0
var _span: float = 1.0
var _unit: String = "s"

const KEYS: Array[String] = ["volume", "lipid_area", "reduced", "catalyst"]
const NAMES: Array[String] = ["Volume (V0)", "Lipid area (A0)", "Reduced volume", "Catalyst (M)"]
const COLOURS: Array[Color] = [Color8(16, 179, 154), Color8(43, 89, 255), Color8(139, 108, 240), Color8(227, 155, 18)]


func _init() -> void:
	custom_minimum_size = Vector2(700, 320)
	for i in KEYS.size():
		_series.append(PackedFloat32Array())


func show_frame(recording: Recording, frame: RecordingFrame) -> void:
	_unit = recording.time_unit
	var count := recording.frame_count()
	while _read < count:
		var f := recording.frame(_read)
		_read += 1
		_times.append(f.time)
		for i in KEYS.size():
			_series[i].append(recording.value(f, KEYS[i]))
		if recording.value(f, "event") == 1.0:
			_divisions.append(f.time)
	_current = frame.index
	_span = maxf(1.0, _times[_times.size() - 1] if not _times.is_empty() else 1.0)
	queue_redraw()


func _draw() -> void:
	draw_style_box(Look.panel(Color(1, 1, 1, 0.93), 16, 0), Rect2(Vector2.ZERO, size))
	draw_string(Look.font("title"), Vector2(20, 30), "Over time (logarithmic scale)", HORIZONTAL_ALIGNMENT_LEFT, -1, 17, Look.INK)
	var x := 20.0
	for i in NAMES.size():
		draw_rect(Rect2(x, 47, 16, 3), COLOURS[i])
		draw_string(Look.font("medium"), Vector2(x + 22, 52), NAMES[i], HORIZONTAL_ALIGNMENT_LEFT, -1, 11, Look.INK)
		x += Look.font("medium").get_string_size(NAMES[i], HORIZONTAL_ALIGNMENT_LEFT, -1, 11).x + 40.0
	var plot := Rect2(54, 70, size.x - 80, size.y - 110)
	draw_rect(plot, Look.LINE, false, 1.0)
	var low := -2.0   # log10
	var high := 3.0
	for mark: int in [-2, -1, 0, 1, 2, 3]:
		var at_y: float = plot.end.y - (mark - low) / (high - low) * plot.size.y
		draw_line(Vector2(plot.position.x, at_y), Vector2(plot.end.x, at_y), Look.LINE, 1.0)
		draw_string(Look.font("mono"), Vector2(10, at_y + 4), Look.significant(pow(10.0, mark), 1), HORIZONTAL_ALIGNMENT_LEFT, -1, 11,
			Look.MUTED)
	for k: int in 5:
		var t: float = _span * k / 4.0
		var at_x := plot.position.x + t / _span * plot.size.x
		var text := Look.time(t, _unit) if _unit == "s" else str(int(round(t)))  # the unit is on the time bar
		draw_string(Look.font("mono"), Vector2(at_x - 40, plot.end.y + 16), text, HORIZONTAL_ALIGNMENT_CENTER, 80, 11,
			Look.MUTED)
	for t: float in _divisions:
		if _times.size() > _current and t > _times[_current]:
			continue
		var at := plot.position.x + t / _span * plot.size.x
		draw_line(Vector2(at, plot.position.y), Vector2(at, plot.end.y), Color(Look.AMBER, 0.35), 1.0)
	for i in _series.size():
		var points := PackedVector2Array()
		for k in mini(_current + 1, _times.size()):
			var v := _series[i][k]
			if is_nan(v) or v <= 0.0:
				continue
			var y := clampf(log(v) / log(10.0), low, high)
			points.append(Vector2(plot.position.x + _times[k] / _span * plot.size.x, plot.end.y - (y - low) / (high - low) * plot.size.y))
		if points.size() > 1:
			draw_polyline(points, COLOURS[i], 2.0, true)

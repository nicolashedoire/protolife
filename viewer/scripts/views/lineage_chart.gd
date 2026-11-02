class_name LineageChart
extends Control
## A lineage of vesicles over time, up to the frame shown, on a logarithmic scale: the descendants of a vesicle of the
## start, the cargo each one holds (against a vesicle of the start's), their reduced volume, and the vesicles born
## empty of the micelles; the shears as bands.

const LOW := 0.005
const HIGH := 500.0

var _times := PackedFloat32Array()
var _series: Array[PackedFloat32Array] = []
var _names: PackedStringArray = []
var _colours: Array[Color] = []
var _read: int = 0
var _current: int = 0
var _span: float = 60.0
var _shears: Array[Vector2] = []


func _init() -> void:
	custom_minimum_size = Vector2(700, 320)


func _prepare(recording: Recording) -> void:
	var sugars := recording.constant("inside_sugars", 0.0) > 0.5
	for item: Array in [["Descendants per starting vesicle", Look.LIVING],
			["Sugars per vesicle (seeded = 1)" if sugars else "Cargo per vesicle (start = 1)", Look.AMBER],
			["Reduced volume", Look.PRIMARY],
			["Born of micelles, per starting vesicle" if sugars else "Vesicles born empty, per starting vesicle",
				Color8(120, 132, 150)]]:
		_names.append(str(item[0]))
		_colours.append(item[1])
		_series.append(PackedFloat32Array())
	if recording.constants.has("cycles") and recording.constant("stirred", 0.0) < 0.5:
		for c in int(recording.constant("cycles")):
			var start := c * recording.constant("period") + recording.constant("shear_at")
			_shears.append(Vector2(start, start + recording.constant("shear_for")))


## The vesicles born of the micelles in a frame (the groups of origin 1 and more), per vesicle of the start.
static func born(recording: Recording, frame: RecordingFrame) -> float:
	var total := 0.0
	for g in range(1, int(recording.constant("records", 1.0))):
		if recording.value(frame, "r%d_present" % g) > 0.5 and recording.value(frame, "r%d_origin" % g) > 0.5:
			total += recording.value(frame, "r%d_count" % g)
	return total


func show_frame(recording: Recording, frame: RecordingFrame) -> void:
	if _series.is_empty():
		_prepare(recording)
	var count := recording.frame_count()
	while _read < count:
		var f := recording.frame(_read)
		_read += 1
		_times.append(f.time)
		_series[0].append(recording.value(f, "r0_count"))
		_series[1].append(recording.value(f, "r0_cargo_share"))
		_series[2].append(recording.value(f, "r0_reduced"))
		_series[3].append(born(recording, f))
	_current = frame.index
	_span = maxf(60.0, _times[_times.size() - 1] if not _times.is_empty() else 60.0)
	queue_redraw()


func _y(plot: Rect2, value: float) -> float:
	var t := (log(clampf(value, LOW, HIGH)) - log(LOW)) / (log(HIGH) - log(LOW))
	return plot.end.y - t * plot.size.y


func _draw() -> void:
	draw_style_box(Look.panel(Color(1, 1, 1, 0.93), 16, 0), Rect2(Vector2.ZERO, size))
	draw_string(Look.font("title"), Vector2(20, 30), "Over time (logarithmic scale)", HORIZONTAL_ALIGNMENT_LEFT, -1, 17,
		Look.INK)
	var x := 20.0
	var y := 52.0
	for i in _names.size():
		var width := Look.font("medium").get_string_size(_names[i], HORIZONTAL_ALIGNMENT_LEFT, -1, 11).x + 34.0
		if x + width > size.x - 20:
			x = 20.0
			y += 16.0
		draw_rect(Rect2(x, y - 5, 16, 3), _colours[i])
		draw_string(Look.font("medium"), Vector2(x + 22, y), _names[i], HORIZONTAL_ALIGNMENT_LEFT, -1, 11, Look.INK)
		x += width
	var plot := Rect2(54, y + 16, size.x - 80, size.y - y - 50)
	for band: Vector2 in _shears:
		var from := plot.position.x + band.x / _span * plot.size.x
		var to := plot.position.x + band.y / _span * plot.size.x
		if from < plot.end.x:
			draw_rect(Rect2(from, plot.position.y, maxf(minf(to, plot.end.x) - from, 2.0), plot.size.y), Color(Look.VIOLET, 0.18))
	draw_rect(plot, Look.LINE, false, 1.0)
	for mark: float in [0.01, 0.1, 1.0, 10.0, 100.0]:
		var at_y := _y(plot, mark)
		draw_line(Vector2(plot.position.x, at_y), Vector2(plot.end.x, at_y), Look.LINE, 1.0)
		draw_string(Look.font("mono"), Vector2(8, at_y + 4), Look.significant(mark, 2), HORIZONTAL_ALIGNMENT_LEFT, -1, 11, Look.MUTED)
	for tick: Dictionary in Look.ticks(_span):
		var at_x := plot.position.x + float(tick["time"]) / _span * plot.size.x
		draw_string(Look.font("mono"), Vector2(at_x - 20, plot.end.y + 16), str(tick["text"]), HORIZONTAL_ALIGNMENT_CENTER, 40, 11,
			Look.MUTED)
	for i in _series.size():
		var points := PackedVector2Array()
		for k in mini(_current + 1, _times.size()):
			var v := _series[i][k]
			if is_nan(v) or v <= 0.0:
				continue
			points.append(Vector2(plot.position.x + _times[k] / _span * plot.size.x, _y(plot, v)))
		if points.size() > 1:
			draw_polyline(points, _colours[i], 2.4, true)

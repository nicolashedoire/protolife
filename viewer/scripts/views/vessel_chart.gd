class_name VesselChart
extends Control
## The vessel over time, up to the frame shown: in a batch, the oil left, the ester the aggregates carried and the
## ester hydrolysed, as shares of the starting ester, with Buhse's 95 min; in a continuous reactor, the oil left, the
## share of the hydrolysed ester the aggregates carried, and the conversion of the ester fed.

var _times := PackedFloat32Array()
var _series: Array[PackedFloat32Array] = [PackedFloat32Array(), PackedFloat32Array(), PackedFloat32Array()]
var _read: int = 0
var _current: int = 0
var _flowing: bool = false
var _span: float = 6600.0


func _init() -> void:
	custom_minimum_size = Vector2(700, 230)


func show_frame(recording: Recording, frame: RecordingFrame) -> void:
	_catch_up(recording)
	_current = frame.index
	_span = maxf(6600.0, frame.time * 1.05)
	queue_redraw()


## Reads the frames not read yet: a vessel's frames are a hundred bytes.
func _catch_up(recording: Recording) -> void:
	var count := recording.frame_count()
	var ester := recording.constant("ester_mol", NAN)
	var first := recording.frame(0)
	var oil0 := recording.value(first, "oil_ml") if first != null else NAN
	while _read < count:
		var f := recording.frame(_read)
		_read += 1
		var hydrolysed := recording.value(f, "hydrolysis")
		var fed := recording.value(f, "inflow_ester")
		_flowing = _flowing or fed > 0.0
		_times.append(f.time)
		_series[0].append(recording.value(f, "oil_ml") / oil0)
		if fed > 0.0:
			var carried := recording.value(f, "transfer") * recording.constant("esters_per_aggregate", 4.5)
			_series[1].append(carried / hydrolysed if hydrolysed > 0.0 else 0.0)
			_series[2].append(hydrolysed / fed)
		else:
			_series[1].append(recording.value(f, "carried_mol") / ester)
			_series[2].append((recording.value(f, "ethanol_mol") - recording.value(first, "ethanol_mol")) / ester)


func _draw() -> void:
	draw_style_box(Look.panel(Color(1, 1, 1, 0.93), 16, 0), Rect2(Vector2.ZERO, size))
	draw_string(Look.font("title"), Vector2(20, 30), "Over time", HORIZONTAL_ALIGNMENT_LEFT, -1, 17, Look.INK)
	var names := ["Oil left", "Share carried by the aggregates" if _flowing else "Ester carried by the aggregates",
		"Conversion of the ester fed" if _flowing else "Ester hydrolysed"]
	var colours := [Look.AMBER, Look.LIVING, Look.PRIMARY]
	var x := 170.0
	for i in 3:
		draw_rect(Rect2(x, 22, 16, 3), colours[i])
		draw_string(Look.font("medium"), Vector2(x + 22, 28), str(names[i]), HORIZONTAL_ALIGNMENT_LEFT, -1, 12, Look.INK)
		x += 32.0 + Look.font("medium").get_string_size(str(names[i]), HORIZONTAL_ALIGNMENT_LEFT, -1, 12).x
	var plot := Rect2(54, 48, size.x - 80, size.y - 82)
	draw_rect(plot, Look.LINE, false, 1.0)
	for p: int in [0, 50, 100]:
		var y := plot.end.y - p / 110.0 * plot.size.y
		draw_line(Vector2(plot.position.x, y), Vector2(plot.end.x, y), Look.LINE, 1.0)
		draw_string(Look.font("mono"), Vector2(8, y + 4), "%d %%" % p, HORIZONTAL_ALIGNMENT_LEFT, -1, 11, Look.MUTED)
	for tick: Dictionary in Look.ticks(_span):
		var at_x := plot.position.x + float(tick["time"]) / _span * plot.size.x
		draw_string(Look.font("mono"), Vector2(at_x - 20, plot.end.y + 18), str(tick["text"]), HORIZONTAL_ALIGNMENT_CENTER, 40,
			11, Look.MUTED)
	if not _flowing:
		var at := plot.position.x + 5700.0 / _span * plot.size.x
		draw_line(Vector2(at, plot.position.y), Vector2(at, plot.end.y), Look.AMBER, 1.5)
		draw_string(Look.font("medium"), Vector2(at + 6, plot.end.y - 10), "Buhse: 95 min", HORIZONTAL_ALIGNMENT_LEFT, -1, 11,
			Look.AMBER_TEXT)
	for s in 3:
		var points := PackedVector2Array()
		for i in mini(_current + 1, _times.size()):
			var share := clampf(_series[s][i], 0.0, 1.1) if not is_nan(_series[s][i]) else 0.0
			points.append(Vector2(plot.position.x + _times[i] / _span * plot.size.x, plot.end.y - share / 1.1 * plot.size.y))
		if points.size() > 1:
			draw_polyline(points, colours[s], 2.2, true)

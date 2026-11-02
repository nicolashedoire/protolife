class_name VesselLoop
extends Control
## The loop of the vessel, drawn over the 3D view: oil → interface → dissolved ester → octanoate, the octanoate and
## the interface's ester gathering into mixed aggregates that carry the ester into the water, and the empty micelles.
## Each stock shows its amount, each arrow its net flow (thicker when larger, reversed when the flow goes back); the
## loop is green. Read from one frame of the recording.

const NODE := Vector2(132, 54)

var _amounts: Dictionary = {}
var _flows: Dictionary = {}
var _flowing: bool = false
var _feed: float = 0.0


func _init() -> void:
	custom_minimum_size = Vector2(700, 320)


func show_frame(recording: Recording, frame: RecordingFrame) -> void:
	for key: String in ["oil_mol", "interface_mol", "ester_mol", "octanoate_mol", "mixed_mol", "micelles_mol"]:
		_amounts[key] = recording.value(frame, key)
	for key: String in ["dispersion", "dissolution", "hydrolysis", "aggregation", "transfer", "micellization"]:
		_flows[key] = recording.value(frame, key)
	_feed = recording.value(frame, "inflow_ester")
	_flowing = _feed > 0.0
	queue_redraw()


func _draw() -> void:
	draw_style_box(Look.panel(Color(1, 1, 1, 0.93), 16, 0), Rect2(Vector2.ZERO, size))
	var font := Look.font("title")
	draw_string(font, Vector2(20, 32), "The loop", HORIZONTAL_ALIGNMENT_LEFT, -1, 19, Look.INK)
	var column := (size.x - 40 - NODE.x) / 3.0
	var top := 96.0
	var spots := {
		"oil_mol": Vector2(20, top), "interface_mol": Vector2(20 + column, top),
		"ester_mol": Vector2(20 + 2 * column, top), "octanoate_mol": Vector2(20 + 3 * column, top),
		"mixed_mol": Vector2(20 + 1.5 * column, top + 128), "micelles_mol": Vector2(20 + 3 * column, top + 128)}
	var names := {"oil_mol": "Oil", "interface_mol": "Interface", "ester_mol": "Dissolved ester",
		"octanoate_mol": "Free octanoate", "mixed_mol": "Mixed aggregates", "micelles_mol": "Empty micelles"}
	var row := top + NODE.y / 2
	_arrow(Vector2(20 + NODE.x, row), Vector2(20 + column, row), "dispersion", "dispersion", "coalescence", false, 0)
	_arrow(Vector2(20 + column + NODE.x, row), Vector2(20 + 2 * column, row), "dissolution", "dissolution",
		"back to the interface", false, 0)
	_arrow(Vector2(20 + 2 * column + NODE.x, row), Vector2(20 + 3 * column, row), "hydrolysis", "hydrolyse", "", true, 0)
	var mixed: Vector2 = spots["mixed_mol"]
	var soap: Vector2 = spots["octanoate_mol"]
	var ester: Vector2 = spots["ester_mol"]
	var face: Vector2 = spots["interface_mol"]
	var corner := Vector2(face.x + 30, mixed.y + NODE.y / 2)
	draw_line(face + Vector2(30, NODE.y), corner, Look.LINE, 1.5)
	draw_line(corner, mixed + Vector2(0, NODE.y / 2), Look.LINE, 1.5)
	draw_string(Look.font("medium"), face + Vector2(38, NODE.y + 30), "4.5 esters", HORIZONTAL_ALIGNMENT_LEFT, -1, 11, Look.MUTED)
	draw_string(Look.font("medium"), face + Vector2(38, NODE.y + 44), "per aggregate", HORIZONTAL_ALIGNMENT_LEFT, -1, 11, Look.MUTED)
	_arrow(soap + Vector2(20, NODE.y), mixed + Vector2(NODE.x, 18), "aggregation", "aggregation", "release", true, 1)
	_arrow(mixed + Vector2(NODE.x / 2 + 26, 0), ester + Vector2(30, NODE.y), "transfer", "transfert", "reformation", true, -1)
	_arrow(soap + Vector2(NODE.x - 24, NODE.y), spots["micelles_mol"] + Vector2(NODE.x - 24, 0), "micellization",
		"micellisation", "dissociation", false, -1)
	for key: String in spots:
		var loop := key in ["ester_mol", "octanoate_mol", "mixed_mol"]
		var box := Look.panel(Look.LIVING_SOFT if loop else Look.SURFACE, 10, 0)
		box.border_color = Look.LIVING if loop else Look.LINE
		box.set_border_width_all(2 if loop else 1)
		draw_style_box(box, Rect2(spots[key], NODE))
		draw_string(Look.font("bold"), spots[key] + Vector2(12, 22), str(names[key]), HORIZONTAL_ALIGNMENT_LEFT, -1, 14, Look.INK)
		draw_string(Look.font("medium"), spots[key] + Vector2(12, 42), Look.moles(float(_amounts.get(key, NAN))),
			HORIZONTAL_ALIGNMENT_LEFT, NODE.x - 16, 13, Look.LIVING_TEXT if loop else Look.MUTED)
	var line := "Closed batch: nothing in, nothing out."
	if _flowing:
		line = "Continuous reactor: ester fed " + Look.moles(_feed, "/min") + ", the mixture carried out at the same rate."
	draw_string(Look.font("medium"), Vector2(20, size.y - 16), line, HORIZONTAL_ALIGNMENT_LEFT, size.x - 40, 13, Look.INK)


## An arrow drawn the way the net flow goes, thicker when larger, with the step's name and rate: above it (side 0),
## on its left (-1) or on its right (1).
func _arrow(from: Vector2, to: Vector2, key: String, forward: String, reverse: String, loop: bool, side: int) -> void:
	var rate := float(_flows.get(key, NAN))
	var size_of := absf(rate) if not is_nan(rate) else 0.0
	var alive := size_of * Look.AVOGADRO >= 1.0
	var colour := (Look.LIVING if loop else Look.PRIMARY) if alive else Look.LINE
	var middle := (from + to) / 2.0
	if rate < 0.0:
		var swap := from
		from = to
		to = swap
	var unit := (to - from).normalized()
	var normal := Vector2(-unit.y, unit.x)
	var tip := to - unit * 3.0
	var base := tip - unit * 13.0
	var width := clampf(log(size_of) / log(10.0) + 9.0, 1.5, 7.0) if alive else 1.0
	draw_line(from, base, colour, width, true)
	draw_colored_polygon(PackedVector2Array([tip, base + normal * 7.0, base - normal * 7.0]), colour)
	var name := reverse if rate < 0.0 and not reverse.is_empty() else forward
	var ink := Look.LIVING_TEXT if loop else Look.MUTED
	var value := Look.moles(size_of, "/min") if alive else "stopped"
	var at := Vector2(middle.x - 80, minf(from.y, to.y) - 46)
	var align := HORIZONTAL_ALIGNMENT_CENTER
	if side < 0:
		at = middle + Vector2(-172, -6)
		align = HORIZONTAL_ALIGNMENT_RIGHT
	elif side > 0:
		at = middle + Vector2(14, 6)
		align = HORIZONTAL_ALIGNMENT_LEFT
	draw_string(Look.font("bold"), at, name, align, 160, 11, ink)
	draw_string(Look.font("medium"), at + Vector2(0, 14), value, align, 160, 11, ink)

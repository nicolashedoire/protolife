class_name Look
extends RefCounted
## The look of the viewer: the colours and fonts of the Genesis lab's screens (where protolife comes from), the
## panels and labels built from them, and numbers written the French way.

const BACKGROUND := Color8(243, 245, 248)
const PANEL := Color8(255, 255, 255)
const LINE := Color8(226, 231, 238)
const INK := Color8(14, 26, 43)
const MUTED := Color8(90, 107, 128)
const PRIMARY := Color8(43, 89, 255)
const PRIMARY_SOFT := Color8(238, 242, 255)
const SURFACE := Color8(247, 249, 252)
const LIVING := Color8(16, 179, 154)
const LIVING_SOFT := Color8(230, 247, 243)
const LIVING_TEXT := Color8(11, 122, 104)
const AMBER := Color8(227, 155, 18)
const AMBER_TEXT := Color8(138, 90, 0)
const VIOLET := Color8(139, 108, 240)
const VIOLET_TEXT := Color8(75, 52, 184)
const DANGER := Color8(194, 54, 47)
const NIGHT := Color8(14, 26, 43)
const NIGHT_MUTED := Color8(154, 168, 188)

const AVOGADRO := 6.02214076e23

static var _fonts: Dictionary = {}


## "body", "medium", "bold", "title" or "mono", read from the TTF files themselves (no import step needed).
static func font(role: String) -> Font:
	if _fonts.has(role):
		return _fonts[role]
	var files := {"body": "Inter-Regular.ttf", "medium": "Inter-Regular.ttf", "bold": "Inter-SemiBold.ttf",
		"title": "Sora-SemiBold.ttf", "mono": "JetBrainsMono-Regular.ttf"}
	var file := FontFile.new()
	file.data = FileAccess.get_file_as_bytes("res://assets/fonts/" + str(files.get(role, "Inter-Regular.ttf")))
	_fonts[role] = file
	return file


static func panel(fill: Color = PANEL, radius: int = 14, padding: int = 18) -> StyleBoxFlat:
	var box := StyleBoxFlat.new()
	box.bg_color = fill
	box.set_corner_radius_all(radius)
	box.set_content_margin_all(padding)
	box.border_color = LINE
	box.set_border_width_all(1 if fill == PANEL else 0)
	return box


static func label(text: String, size: int, color: Color = INK, role: String = "body") -> Label:
	var result := Label.new()
	result.text = text
	result.add_theme_font_override("font", font(role))
	result.add_theme_font_size_override("font_size", size)
	result.add_theme_color_override("font_color", color)
	return result


static func button(text: String, active: bool = false) -> Button:
	var result := Button.new()
	result.text = text
	result.focus_mode = Control.FOCUS_NONE
	result.add_theme_font_override("font", font("bold"))
	result.add_theme_font_size_override("font_size", 15)
	var normal := panel(INK if active else PANEL, 10, 10)
	normal.content_margin_top = 7
	normal.content_margin_bottom = 7
	normal.set_border_width_all(0 if active else 1)
	var hover := normal.duplicate() as StyleBoxFlat
	hover.bg_color = INK if active else PRIMARY_SOFT
	result.add_theme_stylebox_override("normal", normal)
	result.add_theme_stylebox_override("hover", hover)
	result.add_theme_stylebox_override("pressed", hover)
	result.add_theme_color_override("font_color", PANEL if active else INK)
	result.add_theme_color_override("font_hover_color", PANEL if active else PRIMARY)
	result.add_theme_color_override("font_pressed_color", PANEL if active else PRIMARY)
	return result


## A number with `digits` significant digits.
static func significant(value: float, digits: int = 3) -> String:
	if is_nan(value):
		return "—"
	if value == 0.0:
		return "0"
	if absf(value) < 1e-3 or absf(value) >= 1e6:
		return scientific(value, digits)
	var magnitude := floori(log(absf(value)) / log(10.0))
	var decimals := clampi(digits - 1 - magnitude, 0, 6)
	var text := String.num(value, decimals)
	if text.contains("."):
		text = text.rstrip("0").rstrip(".")
	return text


## A number in scientific notation: "1.4·10⁻¹⁴".
static func scientific(value: float, digits: int = 2) -> String:
	if is_nan(value):
		return "—"
	if value == 0.0:
		return "0"
	var exponent := floori(log(absf(value)) / log(10.0))
	var mantissa := value / pow(10.0, exponent)
	var superscript := ""
	for c in str(exponent):
		superscript += str({"-": "⁻", "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷",
			"8": "⁸", "9": "⁹"}[c])
	return String.num(mantissa, digits - 1) + "·10" + superscript


## Ticks for an axis from 0 to `span` seconds: a round step in minutes or hours, at most six ticks.
static func ticks(span: float) -> Array[Dictionary]:
	var steps := [60.0, 300.0, 600.0, 1200.0, 1800.0, 3600.0, 7200.0, 14400.0, 28800.0, 86400.0, 172800.0]
	var step: float = steps[steps.size() - 1]
	for s: float in steps:
		if span / s <= 6.0:
			step = s
			break
	var list: Array[Dictionary] = []
	var t := 0.0
	while t <= span + 1e-6:
		var text := ("%d min" % int(t / 60.0)) if step < 3600.0 else ("%d h" % int(t / 3600.0))
		list.append({"time": t, "text": text})
		t += step
	return list


## An amount with mol, mmol, µmol, ... or "less than one molecule" below one.
static func moles(value: float, suffix: String = "") -> String:
	if is_nan(value):
		return "—"
	var size := absf(value)
	if size * AVOGADRO < 1.0:
		return "less than one molecule"
	var units := [["mol", 1.0], ["mmol", 1e-3], ["µmol", 1e-6], ["nmol", 1e-9], ["pmol", 1e-12], ["fmol", 1e-15],
		["amol", 1e-18], ["zmol", 1e-21]]
	for unit: Array in units:
		var scale: float = unit[1]
		if size >= scale * 0.9995:
			return significant(value / scale) + " " + str(unit[0]) + suffix
	return significant(value / 1e-24) + " ymol" + suffix


## A simulated duration: "1 h 28 min", "42 s", "3 j 4 h".
static func duration(seconds: float) -> String:
	if is_nan(seconds):
		return "—"
	var total := int(round(seconds))
	var days := total / 86400
	var hours := (total % 86400) / 3600
	var minutes := (total % 3600) / 60
	if days > 0:
		return "%d d %d h" % [days, hours]
	if hours > 0:
		return "%d h %02d min" % [hours, minutes]
	if minutes > 0:
		return "%d min %02d s" % [minutes, total % 60]
	return "%d s" % total

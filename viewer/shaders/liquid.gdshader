// The liquid of the biphasic vessel: a sphere cut at the oil's surface, water below the interface, oil above it.
// The levels come from the recorded volumes (scripts/views/vessel_view.gd), in the sphere's own coordinates.
shader_type spatial;
render_mode cull_disabled, depth_draw_always;

uniform float water_top = 0.0;
uniform float oil_top = 0.0;
uniform vec4 water_color : source_color = vec4(0.58, 0.78, 0.93, 0.5);
uniform vec4 oil_color : source_color = vec4(0.90, 0.69, 0.26, 0.85);

varying vec3 local;

void vertex() {
	local = VERTEX;
}

void fragment() {
	if (local.y > oil_top) {
		discard;
	}
	vec4 colour = local.y > water_top ? oil_color : water_color;
	ALBEDO = colour.rgb;
	ALPHA = colour.a;
	ROUGHNESS = 0.18;
	SPECULAR = 0.6;
}

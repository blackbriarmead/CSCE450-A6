#version 120

varying vec3 vPos; // in camera space
varying vec3 vNor; // in camera space
varying vec2 vTex;
uniform vec3 kd;

void main()
{
	
	vec3 terrain_color = kd;
	float incidence = dot(normalize(vPos), normalize(vNor));
	if (incidence < 0.0) {
		incidence = -incidence;
	}

	incidence = 1.0 - incidence;
	//vec3 color = 0.01 * normalize(vPos) + 0.5 * normalize(vNor);
	vec3 surface_color = mix(terrain_color, vec3(0.0, 0.0, 0.0), incidence);
	float distance = length(vPos);
	float a = 5.0;

	float fogAmount = 1-exp(-distance/a);
	vec3 fogColor = vec3(1.0);

	vec3 color = mix(surface_color, fogColor, fogAmount);

	gl_FragColor = vec4(color, 1.0);
}


//#version 120
//
//varying vec3 vPos; // in camera space
//varying vec3 vNor; // in camera space
//varying vec2 vTex;
////uniform vec3 kdFront;
////uniform vec3 kdBack;
//
//void main()
//{
//	vec3 lightPos = vec3(0.0, 0.0, 0.0);
//	vec3 n = normalize(vNor);
//	vec3 l = normalize(lightPos - vPos);
//	vec3 v = -normalize(vPos);
//	vec3 h = normalize(l + v);
//	vec3 kd = vec3(0.0, 0.0, 1.0);
//	float ln = dot(l, n);
//	vec3 diffuse = ln * kd;
//	vec3 color = diffuse;
//	gl_FragColor = vec4(color, 1.0);
//}
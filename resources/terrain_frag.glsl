

#version 120

varying vec3 vPos; // in camera space
varying vec3 vNor; // in camera space
varying vec2 vTex;
uniform vec3 kdFront;
uniform vec3 kdBack;

void main()
{
	vec3 lightPos = vec3(0.0, 0.0, 0.0);
	vec4 ka = vec4(1.0,1.0,1.0,1.0);
	vec3 ks = vec3(1.0,1.0,1.0);
	float s = 10.0f;

	vec3 n = normalize(vNor);
	vec3 l = normalize(lightPos - vPos);
	vec3 v = -normalize(vPos);
	vec3 h = normalize(l + v);
	float blendFactor = vPos.y/100.0f;
	vec4 colorT = mix(
		vec4(1.0,1.0,1.0,1.0),
		vec4(1.0,0.0,0.0,1.0),
		blendFactor
	);
	vec3 colorA = ka.xyz;
	vec3 colorD = max(dot(l, n), 0.0) * colorT.rgb;
	vec3 colorS = pow(max(dot(h, n), 0.0), s) * ks;
	vec3 color = colorA + colorD + colorS;
	//gl_FragColor = vec4(color, colorT.a);
	/*gl_FragColor = mix(
		vec4(1.0,1.0,1.0,1.0),
		vec4(1.0,0.0,0.0,1.0),
		vPos.z/1.0f*/

	gl_FragColor = vec4(1.0f,1.0f,1.0f,1.0f);
}


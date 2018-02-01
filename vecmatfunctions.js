// NOTE: All vectors are being treated as column vectors and
//		 all matrices are treated as being column major
// Vector multiplication was not coded because it's undefined

// Vectors of Size 4
// Vector Addition
function v4_add(v1, v2) {
	var result = [0, 0, 0, 0];
	result[0] = v1[0] + v2[0];
	result[1] = v1[1] + v2[1];
	result[2] = v1[2] + v2[2];
	result[3] = v1[3] + v2[3];
	return result;
}

// Vector Subtraction
function v4_sub(v1, v2) {
	var result = [0, 0, 0, 0];
	result[0] = v1[0] - v2[0];
	result[1] = v1[1] - v2[1];
	result[2] = v1[2] - v2[2];
	result[3] = v1[3] - v2[3];
	return result;
}

// Vector Scalar Multiplication
function v4_scalar_mult(v, s) {
	var result = [0, 0, 0, 0];
	result[0] = v[0] * s;
	result[1] = v[1] * s;
	result[2] = v[2] * s;
	result[3] = v[3] * s;
	return result;
}

// Dot Product
function v4_dot_prod(v1, v2) {
	var result;
	result = (v1[0] * v2[0]) + (v1[1] * v2[1]) + (v1[2] * v2[2]) + (v1[3] * v2[3]);
	return result;
}

// Cross Product
// Since cross product of vec4 is no possible, do cross product using
// x,y,z and set w to 0
function v4_cross_prod(v1, v2) {
	var result = [0, 0, 0, 0];
	result[0] = (v1[1]*v2[2]) - (v1[2]*v2[1]);
	result[1] = -((v1[0]*v2[2]) - (v1[2]*v2[0]));
	result[2] = (v1[0]*v2[1]) - (v1[1]*v2[0]);
	return result;
}

// Magnitude
function v4_length(v1) {
	var magnitude = Math.sqrt(Math.pow(v1[0], 2) + Math.pow(v1[1], 2) + Math.pow(v1[2], 2));
	return magnitude;
}

// Unit Vector
function v4_unit_vec(v1) {
	var magnitude = v4_length(v1);
	var result = v4_scalar_mult(v1, 1 / magnitude);
	return result;
}


// Matrix Functions For 4x4
// Matrix Addition
function m4_add(m1, m2) {
	var result = [m1[0] + m2[0], m1[1] + m2[1], m1[2] + m2[2], m1[3] + m2[3],
				  m1[4] + m2[4], m1[5] + m2[5], m1[6] + m2[6], m1[7] + m2[7],
				  m1[8] + m2[8], m1[9] + m2[9], m1[10] + m2[10], m1[11] + m2[11],
				  m1[12] + m2[12], m1[13] + m2[13], m1[14] + m2[14], m1[15] + m2[15]];
	return result;
}

// Matrix Subtraction
function m4_sub(m1, m2) {
	var result = [m1[0] - m2[0], m1[1] - m2[1], m1[2] - m2[2], m1[3] - m2[3],
				  m1[4] - m2[4], m1[5] - m2[5], m1[6] - m2[6], m1[7] - m2[7],
				  m1[8] - m2[8], m1[9] - m2[9], m1[10] - m2[10], m1[11] - m2[11],
				  m1[12] - m2[12], m1[13] - m2[13], m1[14] - m2[14], m1[15] - m2[15]];
	return result;
}

// Matrix Transpose
function m4_transpose(m) {
	var result = [m[0], m[4], m[8], m[12],
				  m[1], m[5], m[9], m[13],
				  m[2], m[6], m[10], m[14],
				  m[3], m[7], m[11], m[15]];

	return result;
}

// Matrix Multiplication
function m4_mult(m1, m2) {
	var result = m4_identity();
	// Transpose the first matrix because matrix multiplication involves
	// performing dot product on the rows of the left matrix by the columns 
	// of the right matrix. Transposing sets up the matrices in this way.
	// Makes coding the calculation a lot simpler
	m1 = m4_transpose(m1);
	
	var m1x = m1.slice(0, 4);
	var m1y = m1.slice(4, 8);
	var m1z = m1.slice(8, 12);
	var m1w = m1.slice(12, 16);
	
	var m2x = m2.slice(0, 4);
	var m2y = m2.slice(4, 8);
	var m2z = m2.slice(8, 12);
	var m2w = m2.slice(12, 16);
	
	result[0] = v4_dot_prod(m1x, m2x);
	result[1] = v4_dot_prod(m1y, m2x);
	result[2] = v4_dot_prod(m1z, m2x);
	result[3] = v4_dot_prod(m1w, m2x);

	result[4] = v4_dot_prod(m1x, m2y);
	result[5] = v4_dot_prod(m1y, m2y);
	result[6] = v4_dot_prod(m1z, m2y);
	result[7] = v4_dot_prod(m1w, m2y);

	result[8] = v4_dot_prod(m1x, m2z);
	result[9] = v4_dot_prod(m1y, m2z);
	result[10] = v4_dot_prod(m1z, m2z);
	result[11] = v4_dot_prod(m1w, m2z);

	result[12] = v4_dot_prod(m1x, m2w);
	result[13] = v4_dot_prod(m1y, m2w);
	result[14] = v4_dot_prod(m1z, m2w);
	result[15] = v4_dot_prod(m1w, m2w);

	return result;
}

// Matrix Multiplied by Vector
function m4_mult_vec(m, v) {
	var result = [0, 0, 0, 0];
	// Transpose the first matrix because matrix multiplication involves
	// performing dot product on the rows of the left matrix by the column 
	// of the right vector. Transposing sets up the matrices in this way.
	// Makes coding the calculation a lot simpler
	m = m4_transpose(m);
	
	var mx = m.slice(0, 4);
	var my = m.slice(4, 8);
	var mz = m.slice(8, 12);
	var mw = m.slice(12, 16);

	result[0] = v4_dot_prod(mx, v);
	result[1] = v4_dot_prod(my, v);
	result[2] = v4_dot_prod(mz, v);
	result[3] = v4_dot_prod(mw, v);

	return result;
}

// Matrix Scalar Multiplication
function m4_scalar_mult(m, s) {
	var mx = m.slice(0, 4);
	var my = m.slice(4, 8);
	var mz = m.slice(8, 12);
	var mw = m.slice(12, 16);
	
	var v1 = v4_scalar_mult(mx, s);
	var v2 = v4_scalar_mult(my, s);
	var v3 = v4_scalar_mult(mz, s);
	var v4 = v4_scalar_mult(mw, s);
	
	var result = [v1[0], v1[1], v1[2], v1[3],
				  v2[0], v2[1], v2[2], v2[3],
				  v3[0], v3[1], v3[2], v3[3],
				  v4[0], v4[1], v4[2], v4[3]];
	return result;
}

function m4_minor(m) {
	var m11 = (m[5])*(m[10])*(m[15]) + (m[9])*(m[14])*(m[7]) + (m[13])*(m[6])*(m[11]) - (m[7])*(m[10])*(m[13]) - (m[11])*(m[14])*(m[5]) - (m[15])*(m[6])*(m[9]);
	var m12 = (m[4])*(m[10])*(m[15]) + (m[8])*(m[14])*(m[7]) + (m[12])*(m[6])*(m[11]) - (m[7])*(m[10])*(m[12]) - (m[11])*(m[14])*(m[4]) - (m[15])*(m[6])*(m[8]);
	var m13 = (m[4])*(m[9])*(m[15]) + (m[8])*(m[13])*(m[7]) + (m[12])*(m[5])*(m[11]) - (m[7])*(m[9])*(m[12]) - (m[11])*(m[13])*(m[4]) - (m[15])*(m[5])*(m[8]);
	var m14 = (m[4])*(m[9])*(m[14]) + (m[8])*(m[13])*(m[6]) + (m[12])*(m[5])*(m[10]) - (m[6])*(m[9])*(m[12]) - (m[10])*(m[13])*(m[4]) - (m[14])*(m[5])*(m[8]);

	var m21 = (m[1])*(m[10])*(m[15]) + (m[9])*(m[14])*(m[3]) + (m[13])*(m[2])*(m[11]) - (m[3])*(m[10])*(m[13]) - (m[11])*(m[14])*(m[1]) - (m[15])*(m[2])*(m[9]);
	var m22 = (m[0])*(m[10])*(m[15]) + (m[8])*(m[14])*(m[3]) + (m[12])*(m[2])*(m[11]) - (m[3])*(m[10])*(m[12]) - (m[11])*(m[14])*(m[0]) - (m[15])*(m[2])*(m[8]);
	var m23 = (m[0])*(m[9])*(m[15]) + (m[8])*(m[13])*(m[3]) + (m[12])*(m[1])*(m[11]) - (m[3])*(m[9])*(m[12]) - (m[11])*(m[13])*(m[0]) - (m[15])*(m[1])*(m[8]);
	var m24 = (m[0])*(m[9])*(m[14]) + (m[8])*(m[13])*(m[2]) + (m[12])*(m[1])*(m[10]) - (m[2])*(m[9])*(m[12]) - (m[10])*(m[13])*(m[0]) - (m[14])*(m[1])*(m[8]);

	var m31 = (m[1])*(m[6])*(m[15]) + (m[5])*(m[14])*(m[3]) + (m[13])*(m[2])*(m[7]) - (m[3])*(m[6])*(m[13]) - (m[7])*(m[14])*(m[1]) - (m[15])*(m[2])*(m[5]);
	var m32 = (m[0])*(m[6])*(m[15]) + (m[4])*(m[14])*(m[3]) + (m[12])*(m[2])*(m[7]) - (m[3])*(m[6])*(m[12]) - (m[7])*(m[14])*(m[0]) - (m[15])*(m[2])*(m[4]);
	var m33 = (m[0])*(m[5])*(m[15]) + (m[4])*(m[13])*(m[3]) + (m[12])*(m[1])*(m[7]) - (m[3])*(m[5])*(m[12]) - (m[7])*(m[13])*(m[0]) - (m[15])*(m[1])*(m[4]);
	var m34 = (m[0])*(m[5])*(m[14]) + (m[4])*(m[13])*(m[2]) + (m[12])*(m[1])*(m[6]) - (m[2])*(m[5])*(m[12]) - (m[6])*(m[13])*(m[0]) - (m[14])*(m[1])*(m[4]);

	var m41 = (m[1])*(m[6])*(m[11]) + (m[5])*(m[10])*(m[3]) + (m[9])*(m[2])*(m[7]) - (m[3])*(m[6])*(m[9]) - (m[7])*(m[10])*(m[1]) - (m[11])*(m[2])*(m[5]);
	var m42 = (m[0])*(m[6])*(m[11]) + (m[4])*(m[10])*(m[3]) + (m[8])*(m[2])*(m[7]) - (m[3])*(m[6])*(m[8]) - (m[7])*(m[10])*(m[0]) - (m[11])*(m[2])*(m[4]);
	var m43 = (m[0])*(m[5])*(m[11]) + (m[4])*(m[9])*(m[3]) + (m[8])*(m[1])*(m[7]) - (m[3])*(m[5])*(m[8]) - (m[7])*(m[9])*(m[0]) - (m[11])*(m[1])*(m[4]);
	var m44 = (m[0])*(m[5])*(m[10]) + (m[4])*(m[9])*(m[2]) + (m[8])*(m[1])*(m[6]) - (m[2])*(m[5])*(m[8]) - (m[6])*(m[9])*(m[0]) - (m[10])*(m[1])*(m[4]);
	
	var minor = [m11, m12, m13, m14,
				 m21, m22, m23, m24,
				 m31, m32, m33, m34,
				 m41, m42, m43, m44];
	return minor;
}

// Find the cofactor of a matrix
function m4_cofactor(m) {
	var cofactor = [ m[0], -m[1], m[2], -m[3],
					-m[4], m[5], -m[6], m[7],
					m[8], -m[9], m[10], -m[11],
					-m[12], m[13], -m[14], m[15]];

	return cofactor;
}

// Calculate determinant of matrix
// If determinant returns 0 then an inverse doesn't exist
function m4_determinant(m) {
	var m11 = (m[5])*(m[10])*(m[15]) + (m[9])*(m[14])*(m[7]) + (m[13])*(m[6])*(m[11]) - (m[7])*(m[10])*(m[13]) - (m[11])*(m[14])*(m[5]) - (m[15])*(m[6])*(m[9]);
	var m21 = (m[1])*(m[10])*(m[15]) + (m[9])*(m[14])*(m[3]) + (m[13])*(m[2])*(m[11]) - (m[3])*(m[10])*(m[13]) - (m[11])*(m[14])*(m[1]) - (m[15])*(m[2])*(m[9]);
	var m31 = (m[1])*(m[6])*(m[15]) + (m[5])*(m[14])*(m[3]) + (m[13])*(m[2])*(m[7]) - (m[3])*(m[6])*(m[13]) - (m[7])*(m[14])*(m[1]) - (m[15])*(m[2])*(m[5]);
	var m41 = (m[1])*(m[6])*(m[11]) + (m[5])*(m[10])*(m[3]) + (m[9])*(m[2])*(m[7]) - (m[3])*(m[6])*(m[9]) - (m[7])*(m[10])*(m[1]) - (m[11])*(m[2])*(m[5]);

	var determinant = m[0]*m11 - m[4]*m21 + m[8]*m31 - m[12]*m41;

	return determinant;
}

// Matrix Inverse
function m4_inverse(m) {
	var minor = m4_minor(m);
	var cofactor = m4_cofactor(minor);
	var transpose = m4_transpose(cofactor);
	var determinant = m4_determinant(m);
	var inverse = m4_scalar_mult(transpose, 1 / determinant);
	return inverse;
}

// Identity Matrix
function m4_identity() {
	var result = [1.0, 0.0, 0.0, 0.0,
				  0.0, 1.0, 0.0, 0.0,
				  0.0, 0.0, 1.0, 0.0,
				  0.0, 0.0, 0.0, 1.0];
	return result;
}
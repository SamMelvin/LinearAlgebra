//general mathmatic functions
function sqrt(number) {
  let sqrt = Math.sqrt(number);
  return sqrt;
}

function pow(number, power) {
  let p0wer = Math.pow(number, power);
  return p0wer;
}

//Section 1: Operations on one matrix

function matrixMult(matrix, scalar) {
  //loop through each row of the matrix, multiply each entry by the inputted scalar value.
  for (let i = 0; i < matrix.length; i++) {
    for (let j = 0; j < matrix[0].length; j++) {
      matrix[i][j] = matrix[i][j] * scalar;
    }
  }
  return matrix;
}

function addVectors(scalar1, v1, scalar2, v2) {
  //define a variable for the new, empty vector
  let newVector = [];

  for (let i = 0; i < v1.length; i++) {
    for (let i = 0; i < v1.length; i++) {
      newVector[i] = v1[i] * scalar1 + v2[i] * scalar2;
    }
  }

  return newVector;
}

function subtractVectors(scalar1, v1, scalar2, v2) {
  //define a variable for the new, empty vector
  let newVector = [];

  for (let i = 0; i < v1.length; i++) {
    for (let i = 0; i < v1.length; i++) {
      newVector[i] = v1[i] * scalar1 - v2[i] * scalar2;
    }
  }

  return newVector;
}

function vectorMult(vector, scalar) {
  for (let i = 0; i < vector.length; i++) {
    vector[i] = vector[i] * scalar;
  }
  return vector;
}

function vectorDiv(vector, scalar) {
  let newVector = [];
  for (let i = 0; i < vector.length; i++) {
    newVector[i] = vector[i] / scalar;
  }
  return newVector;
}

//define a function to perform a row operation that zeroes out the pivot entry of Row1 by subratcting a scaled multiple of Row2.
function rowOpSubtract(r1, r2) {
  let r2PivIndex = r2.indexOf(1);
  let coeffRatio = r1[r2PivIndex] / r2[r2PivIndex];
  let newVector = subtractVectors(1, r1, coeffRatio, r2);

  return newVector;
}

//Define a function to divide a row by a scalar, to reduce the pivot entry to 1.
function rowPivTo1(row) {
  //Find the index value of the first nonzero value in the row. Then index this element from the row to find the pivot value.

  function nonZero(value) {
    return value !== 0;
  }

  let pivIndex = row.findIndex(nonZero);
  let pivEntry = row[pivIndex];

  //divide the row by the value of the pivot entry, to reduce the value of the pivot entry to 1.
  let newRow = vectorDiv(row, pivEntry);

  //return the new row, divided by a scalar to reduce the pivot entry to 1.
  return newRow;
}

//needs work
function rref(matrix) {
  //define functions to use
  function rowOpSubtract(r1, r2) {
    let r2PivIndex = r2.indexOf(1);
    let coeffRatio = r1[r2PivIndex] / r2[r2PivIndex];
    let newVector = subtractVectors(1, r1, coeffRatio, r2);

    return newVector;
  }
  //make sure first entry in the first column in nonzero. If not, loop through the matrix until you find a column with a nonzero first entry, and swap the 2 rows. Then, if the new first entry of row 1 is not 1, divide row 1 by a scalar until this is the case.

  if (matrix[0][0] === 0) {
    for (let i = 1; i < matrix.length; i++) {
      if (matrix[0][0] !== 0) {
        break;
      } else if (matrix[i][0] !== 0) {
        let fRow = matrix[0];
        matrix[0] = matrix[i];
        matrix[i] = fRow;
      }
    }
  }

  if (matrix[0][0] !== 1) {
    matrix[0] = matrix[0] / matrix[0][0];
  }

  //create a loop for the length of the matrix
  for (let i = 0; i < matrix.length; i++) {
    //create a loop for the length of each row
    for (let j = 0; j < matrix[0].length; j++) return matrix;
  }
}

//S4: Dot products and Cross Products
function vectorDotProduct(v1, v2) {
  let dotProduct = 0;
  if (v1.length === v2.length) {
    for (let i = 0; i < v1.length; i++) {
      dotProduct = dotProduct + v1[i] * v2[i];
    }
  }

  return dotProduct;
}

function matrixMult(m1, m2) {
  //define the dimensions of the new matrix
  let newMatrix = [];
  for (let i = 0; i < m1.length; i++) {
    newMatrix[i] = [];
  }

  //create a loop for each row of the new matrix
  for (let i = 0; i < m1.length; i++) {
    //create a loop for each entry in each row
    for (let j = 0; j < m2[0].length; j++) {
      //define the row vector from the first matrix
      let m1RV = [];
      for (let k = 0; k < m1[0].length; k++) {
        m1RV[k] = m1[i][k];
      }
      //define the column vector from the second matrix, to multiply the row from the first matrix by
      let m2CV = [];
      for (let k = 0; k < m2.length; k++) {
        m2CV[k] = m2[k][j];
      }

      newMatrix[i][j] = vectorDotProduct(m1RV, m2CV);
    }
  }
  return newMatrix;
}

let m1 = [
  [1, 2, 1],
  [1, 4, 7],
];
let m2 = [
  [1, 2, 3],
  [1, 2, 1],
  [3, 1, 3],
];

//SECTION 5: MATRIX-VECTOR PRODUCTS

//define a function to calculate the length of a given vector.
function vectorLength(vector) {
  //loop through the vector, squaring each value in the vector.
  let sqrtVector = [];
  for (let i = 0; i < vector.length; i++) {
    sqrtVector[i] = pow(vector[i], 2);
  }

  //loop through vector of squared entries, and take the sum of all the entries.
  let sqrtSum = 0;
  for (let i = 0; i < vector.length; i++) {
    sqrtSum = sqrtSum + sqrtVector[i];
  }

  let vectorLength = sqrt(sqrtSum);

  return vectorLength;
}

//SECTION 6: TRANSFORMATIONS
function vectorTransformation(matrix, vector) {
  //define a variable for the new vector
  let newVector = [];

  //loop through the transformation matrix, calculate the dot product for each term in the new vector
  for (let i = 0; i < matrix.length; i++) {
    newVector[i] = vectorDotProduct(matrix[i], vector);
  }

  return newVector;
}

//SECTION 7: INVERSES

//SECTION 8: DETERMINANTS
function crossProduct2x2(matrix) {
  let dotProduct = matrix[0][0] * matrix[1][1] - matrix[1][0] * matrix[0][1];
  return dotProduct;
}

function crossProduct3x3(matrix) {
  let crossProduct =
    matrix[0][0] * (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) -
    matrix[0][1] * (matrix[1][0] * matrix[2][2] - matrix[1][2] * matrix[2][0]) +
    matrix[0][2] * (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]);

  return crossProduct;
}

//define a function to put a square matrix into Upper-Triangular Format.

function upperTriangularFormat(matrix) {
  //make sure that the first entry in row 1 is a nonzero value. If not, swap 2 rows so that this is the case

  //create a loop for each row
  for (let i = 0; i < matrix.length; i++) {}
}

//SECTION 9: TRANSPOSES

//define a function to return the transpose of a matrix.
function transpose(matrix) {
  //define a varibale for the new transpose matrix
  let transpose = [];

  //define dimensions of the transpose matrix
  for (let i = 0; i < matrix[0].length; i++) {
    transpose[i] = [];
  }

  //Define reach row of the transpose, by looping through the matix.

  for (let i = 0; i < transpose.length; i++) {
    for (let j = 0; j < matrix.length; j++) {
      transpose[i][j] = matrix[j][i];
    }
  }

  return transpose;
}

//SECTION 12: ORTHONORMAL BASES AND GRAM-SCHMIDT

//define a function to determine if a given basis is orthonormal or not (returns True or False).
function isOrthonormal(basis) {
  //verify that every vector in the span of the basis has a length of 1
  let allVectorsLength1 = true;
  for (let i = 0; i < basis.length; i++) {
    if (vectorLength(basis[i]) !== 1) {
      allVectorsLength1 = false;
      break;
    }
  }

  //if the length of all of the vectors in the basis are equal to one, confirm that all vector of the new basis are orthogonal to one another.
  let dotResults = [];
  for (let i = 0; i < basis.length; i++) {
    dotResults[i] = [];
  }

  if (allVectorsLength1 === true) {
    let allVectorsOrthogonal = true;
    for (let i = 0; i < basis.length; i++) {
      for (let j = 0; j < basis.length; j++) {
        dotResults[i][j] = 1;
        //dotResults[i][j] = vectorDotProduct(basis[i], basis[j]);
      }

      //dotResults[i][i] = 0;
    }
  }
  dotResults[1] = [2, 2];

  //figure out why the nested loops aren't adding values tonthe 'dotResults' nested array.
  return dotResults;
}

//define a function to calculate the projection of a vector x onto a subspace defined by an orthonormal basis.
function orthonormalBasisProjVX(basisMatrix, vector) {
  //define a variable for the vector projection
  let newVector = [];

  //calculate the new vector with the formula for the projection of a vector X onto a subspace defined by an orthonormal basis
  let ATransposeA = matrixMult(basisMatrix, transpose(basisMatrix));
  let projVX = vectorTransformation(ATransposeA, vector);

  //return the vector projection onto the subspace V

  return projVX;
}

let V = [
  [0, 1 / sqrt(2)],
  [0, 1 / sqrt(2)],
  [1, 0],
];

let W = [
  [0, 0, 1],
  [1 / sqrt(2), 1 / sqrt(2), 0],
];

console.log(rowPivTo1([0, 3, 7, 5]));

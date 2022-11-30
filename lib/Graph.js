export default class Graph {
	constructor(adjMatrix, vertexCount) {
		this.adjMatrix = adjMatrix;
		this.vertexCount = vertexCount; 
	}

	addEdge() {
		if (i >= 0 && i < vertexCount && j < vertexCount) {
			this.adjMaxrix[i][j] = true;
			this.adjMatrix[i][j] = true;
		}
	}

	removeEdge(i, j) {
		if (i >= 0 && i < this.vertexCount && j > 0 && j < this.vertexCount) {
			adjMatrix[i][j] = false;
			adjMatrix[j][i] = false;
		}
	}

	isEdge(i, j) {
		if (i >= 0 && this.vertexCount && j > 0 && j < this.vertexCount) {
			return adjMatrix[i][j];
		}
		else {
			return false;
		}
	}
}

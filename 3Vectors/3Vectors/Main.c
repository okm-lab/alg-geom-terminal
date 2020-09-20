#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
//#define VECTORS_COUNT 3
struct vect {
	int x;
	int y;
	int z;
};

//Функция принимает на вход 2 вектора и возвращает их векторное произведение
struct vect getVectorProduct(struct vect a, struct vect b) {
	struct vect v = { a.y * b.z - a.z* b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x };
	return v;
}

//Функция принимает на вход 2 вектора и возвращает их скалярное произведение
int getScalarProduct(struct vect a, struct vect b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

//Проверка коллинеарности векторов
int checkCollinear(struct vect a, struct vect b) {
	struct vect v = getVectorProduct(a, b);
	return v.x == 0 && v.y == 0 && v.z == 0 ? 1 : 0;
}


int main() {
	//Массив, в котором будут храниться вектора (<x,y,z>)
	struct vect vectors[100];
	//Массив, в котором будут храниться пары векторов, компланарных по отношению к первому
	int coplanar[100][100][2];
	//Массив, в котором будут храниться все вектора, компланарные по отношению к первому
	int collinear[100][100];
	int x1, x2, y1, y2, z1, z2, VECTORS_COUNT;

	printf("Enter number of vectors: ");
	scanf("%d", &VECTORS_COUNT);
	system("cls");

	for (int i = 0; i < VECTORS_COUNT; i++) {
		for (int j = 0; j < VECTORS_COUNT; j++) {
			collinear[i][j] = -1;
			coplanar[i][j][0] = -1;
			coplanar[i][j][1] = -1;
		}
	}
	
	for (int i = 0; i < VECTORS_COUNT; i++) {
		printf("Enter <x1, y1, z1> and <x2, y2, z2> for %d vector: ", i + 1);
		scanf("%d%d%d%d%d%d", &x1, &y1, &z1, &x2, &y2, &z2);
		printf("\n");
		struct vect v = { x2-x1, y2-y1, z2-z1 };
		vectors[i] = v;
		printf("\n");
	}

	for (int i = 0; i < VECTORS_COUNT; i++) {
		//kcomp - индекс следующей пары векторов, компланарных с Vi
		int kcomp = 0;
		//kcol - индекс следующего вектора, коллинеарного Vi
		int kcol = 0;

		for (int j = i + 1; j < VECTORS_COUNT; j++) {
			for (int k = j + 1; k < VECTORS_COUNT; k++) {
				if (getScalarProduct(vectors[i], getVectorProduct(vectors[j], vectors[k])) == 0) {
					coplanar[i][kcomp][0] = j;
					coplanar[i][kcomp][1] = k;
					kcomp++;		
				}
			}
			if (checkCollinear(vectors[i], vectors[j]) == 1) {
				//Проверка того, что этот элемент ещё не встречался в списке коллинеарных векторов
				int chk = 0;
				for (int x = 0; x < j; x++) {
					for (int z = 0; z < VECTORS_COUNT; z++) {
						if (collinear[x][z] == -1)
							break;
						if (collinear[x][z] == j) {
							chk = 1;
							break;
						}
					}
					if (chk == 1)
						break;
				}
				if (chk == 0) {
					collinear[i][kcol] = j;
					kcol++;
				}
			}
		}
	}
	
	for (int i = 0; i < VECTORS_COUNT; i++) {
		if (coplanar[i][0][0] == -1)
			continue;
		for (int j = 0; j < VECTORS_COUNT; j++) {
			if (coplanar[i][j][0] == -1)
				break;
			printf("%d is coplanar with %d and %d\n", i + 1, coplanar[i][j][0]+1, coplanar[i][j][1]+1);
		}
	}

	for (int i = 0; i < VECTORS_COUNT; i++) {
		if (collinear[i][0] == -1)
			continue;
		printf("%d is collinear with: ", i+1);
		int k = 0;
		for (; k < VECTORS_COUNT; k++) {
			if (collinear[i][k] == -1)
				break;
			printf("%d", collinear[i][k] + 1);
			if (k + 1 < VECTORS_COUNT && collinear[i][k + 1] != -1)
				printf(", ");
			else
				printf("\n");
		}
	}
}
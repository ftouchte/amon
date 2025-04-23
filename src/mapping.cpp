/***********************************************
 * mapping.cpp
 *
 * Convert a csv file to c++ routine that 
 * realises the mapping
 *
 * @author Felix Touchte Codjo
 * @date April 21, 2025
 * ********************************************/

#include <cstdio>

int main (int argc, const char* argv[]) {
	const char *filename = "../ressources/ahdc_mapping.csv";
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		perror("Error opening file\n");
		return 1;
	}
	char* hv_group = nullptr;
	int superlayer, layer, wire;
	int extra1, extra2, extra3, extra4, extra5;
	int rows = 0;
	if (fscanf(file, "%s %s %s %s %s %s %s %s %s\n", s1, s2, s3, s4, s5, s6, s7, s8, s9)) { 
		printf();
	}; 
	if (fscanf(file, "%d,%d,%d,%s,%d,%d,%d,%d,%d\n", &superlayer, &layer, &wire, hv_group, &extra1, &extra2, &extra3, &extra4, &extra5) == 9) {
		rows++;
		printf("%d %d %d %s %d %d %d %d %d\n", superlayer, layer, wire, hv_group, extra1, extra2, extra3, extra4, extra5);
	}
	if (fscanf(file, "%d,%d,%d,%s,%d,%d,%d,%d,%d\n", &superlayer, &layer, &wire, hv_group, &extra1, &extra2, &extra3, &extra4, &extra5) == 9) {
		rows++;
		printf("%d %d %d %s %d %d %d %d %d\n", superlayer, layer, wire, hv_group, extra1, extra2, extra3, extra4, extra5);
	}
	while (fscanf(file, "%d,%d,%d,%s,%d,%d,%d,%d,%d\n", &superlayer, &layer, &wire, hv_group, &extra1, &extra2, &extra3, &extra4, &extra5) == 9) {
		rows++;
		printf("%d %d %d %s %d %d %d %d %d\n", superlayer, layer, wire, hv_group, extra1, extra2, extra3, extra4, extra5);
	}
	printf("rows : %d\n", rows);
	return 0;	
}

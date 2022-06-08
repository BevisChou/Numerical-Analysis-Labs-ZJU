OUT = main

lab1: src/lab1.c
	gcc $^ -o ${OUT}

lab6: src/lab6.c
	gcc $^ -o ${OUT}

lab7: src/lab7.c
	gcc $^ -o ${OUT}

lab8: src/lab8.c
	gcc $^ -o ${OUT}

clean:
	$(shell rm ${OUT})
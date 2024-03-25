#create executable for one-root davidson
executable: davidson_template.o driver.o 
	g++ -o executable davidson_template.o driver.o -larmadillo
	echo "The executable was made"

driver.o: driver.cpp
	g++ -o driver.o -c driver.cpp -larmadillo
	echo "The driver object file was made"

davidson_template.o: davidson_template.C
	g++ -o davidson_template.o -c davidson_template.C -larmadillo
	echo "The davidson_template object file was made"

#create debugging executable of the one-root davidson
executable_g: davidson_template_g.o driver_g.o 
	g++ -g -o executable_g davidson_template.o driver.o -larmadillo
	echo "The executable_g was made"

driver_g.o: driver.cpp
	g++ -g -o driver.o -c driver.cpp -larmadillo
	echo "The driver object file was made"

davidson_template_g.o: davidson_template.C
	g++ -g -o davidson_template.o -c davidson_template.C -larmadillo
	echo "The davidson_template object file was made"

#create executable for multi-root davidson
executable_mr: davidson_multiple_roots.o driver_mr.o 
	g++ -o executable_mr davidson_multiple_roots.o driver_mr.o -larmadillo
	echo "The executable was made"

driver_mr.o: driver.cpp
	g++ -o driver_mr.o -c driver.cpp -larmadillo
	echo "The driver object file was made"

davidson_multiple_roots.o: davidson_multiple_roots.C
	g++ -o davidson_multiple_roots.o -c davidson_multiple_roots.C -larmadillo
	echo "The davidson_multiple_roots object file was made"

#create debugging executable for multi-root davidson
executable_mr_g: davidson_multiple_roots_g.o driver_mr_g.o 
	g++ -g -o executable_mr_g davidson_multiple_roots_g.o driver_mr_g.o -larmadillo
	echo "The executable_g was made"

driver_mr_g.o: driver.cpp
	g++ -g -o driver_mr_g.o -c driver.cpp -larmadillo
	echo "The driver object file was made"

davidson_multiple_roots_g.o: davidson_multiple_roots.C
	g++ -g -o davidson_multiple_roots_g.o -c davidson_multiple_roots.C -larmadillo
	echo "The davidson_multiple_roots object file was made"

#create executable for davidson class
executable_class: davidson_class.o driver_class.o 
	g++ -o executable_class davidson_class.o driver_class.o -larmadillo
	echo "The executable was made"

driver_class.o: driver_for_dav_class.cpp
	g++ -o driver_class.o -c driver_for_dav_class.cpp -larmadillo
	echo "The driver object file was made"

davidson_class.o: davidson_class.C
	g++ -o davidson_class.o -c davidson_class.C -larmadillo
	echo "The davidson_class object file was made"

#create debugging executable for davidson class
executable_class_g: davidson_class_g.o driver_class_g.o 
	g++ -g -o executable_class_g davidson_class_g.o driver_class_g.o -larmadillo
	echo "The executable_g was made"

driver_class_g.o: driver_for_dav_class.cpp
	g++ -g -o driver_class_g.o -c driver_for_dav_class.cpp -larmadillo
	echo "The driver object file was made"

davidson_class_g.o: davidson_class.C
	g++ -g -o davidson_class_g.o -c davidson_class.C -larmadillo
	echo "The davidson_class object file was made"


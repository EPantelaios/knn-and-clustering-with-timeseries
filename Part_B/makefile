CC=g++ -std=c++14
CFLAGS=-c -O2
SOURCES_MAIN=Common/main.cpp
SOURCES_CURVE=Curve/Clustering_Frechet.cpp Curve/Discrete.cpp
SOURCES_VECTOR=Vector/Clustering.cpp Vector/Hypercube.cpp Vector/LSH.cpp
SOURCES_COMMON=Common/G_hash.cpp Common/H_hash.cpp Common/Point_frechet.cpp Common/utils.cpp
SOURCES_TESTING=testing.cpp

OBJECTS_CURVE=$(SOURCES_CURVE:.cpp=.o)
OBJECTS_COMMON=$(SOURCES_COMMON:.cpp=.o)
OBJECTS_VECTOR=$(SOURCES_VECTOR:.cpp=.o)
OBJECTS_MAIN=$(SOURCES_MAIN:.cpp=.o)
OBJECTS_TESTING=$(SOURCES_TESTING:.cpp=.o)

EXEC_CURVE=cluster
EXEC_TEST=testing

all: Cluster Cluster_testing
Cluster: $(EXEC_CURVE)
Cluster_testing: $(EXEC_TEST)


$(EXEC_CURVE): $(OBJECTS_MAIN) $(OBJECTS_COMMON) $(OBJECTS_CURVE) $(OBJECTS_VECTOR)
	$(CC) $^ -o $@

$(EXEC_TEST): $(OBJECTS_COMMON) $(OBJECTS_CURVE) $(OBJECTS_VECTOR) $(OBJECTS_TESTING)
	$(CC) $^ -o $@


.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f $(EXEC_CURVE) $(EXEC_TEST) $(OBJECTS_MAIN) $(OBJECTS_CURVE) $(OBJECTS_VECTOR) $(OBJECTS_COMMON) $(OBJECTS_TESTING)

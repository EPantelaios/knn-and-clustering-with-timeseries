CC=g++ -std=c++14
CFLAGS=-c -O2

SOURCES_MAIN=Common/main.cpp
SOURCES_COMMON=Common/utils.cpp Common/G_hash.cpp Common/H_hash.cpp
SOURCES_LSH=LSH/Point_lsh.cpp LSH/LSH.cpp
SOURCES_HYPERCUBE=Hypercube/Hypercube.cpp Hypercube/Point_hc.cpp
SOURCES_DISCRETE=Discrete_Frechet/Curve_discrete.cpp Discrete_Frechet/Discrete.cpp Discrete_Frechet/Point_discrete.cpp
SOURCES_CONTINUOUS=Continuous_Frechet/Continuous.cpp Continuous_Frechet/Curve_cont.cpp Continuous_Frechet/Point_cont.cpp Continuous_Frechet/Fred/src/config.cpp Continuous_Frechet/Fred/src/curve.cpp Continuous_Frechet/Fred/src/frechet.cpp Continuous_Frechet/Fred/src/interval.cpp Continuous_Frechet/Fred/src/point.cpp Continuous_Frechet/Fred/src/simplification.cpp
SOURCES_KNN_TESTING=Unit_Testing/unit_test.cpp

OBJECTS_MAIN=$(SOURCES_MAIN:.cpp=.o)
OBJECTS_COMMON=$(SOURCES_COMMON:.cpp=.o)
OBJECTS_LSH=$(SOURCES_LSH:.cpp=.o)
OBJECTS_HYPERCUBE=$(SOURCES_HYPERCUBE:.cpp=.o)
OBJECTS_DISCRETE=$(SOURCES_DISCRETE:.cpp=.o)
OBJECTS_CONTINUOUS=$(SOURCES_CONTINUOUS:.cpp=.o)
OBJECTS_KNN_TESTING=$(SOURCES_KNN_TESTING:.cpp=.o)

EXEC_KNN=search
EXEC_KNN_TESTING=search_testing


all: knn knn_unit_testing
knn: $(EXEC_KNN)
knn_unit_testing: $(EXEC_KNN_TESTING)


$(EXEC_KNN): $(OBJECTS_MAIN) $(OBJECTS_COMMON) $(OBJECTS_LSH) $(OBJECTS_HYPERCUBE) $(OBJECTS_DISCRETE) $(OBJECTS_CONTINUOUS)
	$(CC) $^ -o $@

$(EXEC_KNN_TESTING): $(OBJECTS_COMMON) $(OBJECTS_LSH) $(OBJECTS_HYPERCUBE) $(OBJECTS_DISCRETE) $(OBJECTS_CONTINUOUS) $(OBJECTS_KNN_TESTING)
	$(CC) $^ -o $@ 


.cpp.o:
	$(CC) $(CFLAGS) $< -o $@


clean:
	rm -f $(EXEC_KNN) $(EXEC_KNN_TESTING) $(OBJECTS_MAIN) $(OBJECTS_COMMON) $(OBJECTS_LSH) $(OBJECTS_HYPERCUBE) $(OBJECTS_DISCRETE) $(OBJECTS_CONTINUOUS) $(OBJECTS_KNN_TESTING)

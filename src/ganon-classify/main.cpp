//CODIGO ORIGINAL
/*#include <ganon-classify/CommandLineParser.hpp>
#include <ganon-classify/GanonClassify.hpp>

#include <cstdlib>
#include <utility>

int main( int argc, char** argv )
{
    if ( auto config = GanonClassify::CommandLineParser::parse( argc, argv ); config.has_value() )
    {
        return GanonClassify::run( std::move( config.value() ) ) ? EXIT_SUCCESS : EXIT_FAILURE;
    }
    else
    {
        return argc == 1 ? EXIT_FAILURE : EXIT_SUCCESS;
    }
}
*/
//FIN CODIGO ORIGINAL


//// CODIGO MODIFICADO PARA USAR MPI, basicmanete la estrategia va a ser dividir el archivo de entrada en partes iguales entre los ranks, y cada rank va a procesar su parte del archivo. Al final, el rank 0 va a fusionar los resultados de todos los ranks. Cada rank usara ganon-classify para procesar su parte del archivo de entrada, y los resultados se guardaran en archivos separados por rank. Este codigo corresponde a una primera estrategia y aproximacion de paralelizacion usando MPI y dejando la implementacion de hilos de C++ base pero usando de momento solo 1 hilo por rank para evitar interferencias entre los procesos MPI.

#include <ganon-classify/CommandLineParser.hpp>
#include <ganon-classify/GanonClassify.hpp>
#include <cstdlib>
#include <utility>

#include <mpi.h>
#include <iostream>
#include <fstream>
#include <filesystem>

namespace fs = std::filesystem;

void merge_files(const std::string& output_prefix_base, const std::string& extension, int num_ranks) //funcion para fusionar los archivos de salida de cada rank a uno solo
{
	//abrimos el archivo de salida (outfile) para escribir, se le llama con el output_prefix_base y la extension correspondiente, por ejemplo ".rep", ".lca", etc. 
	//se le dice que es de salida (std::ios::out) y que es binario (std::ios::binary) para evitar problemas de formato al fusionar los archivos de salida de cada rank
	std::ofstream outfile(output_prefix_base + extension, std::ios::out | std::ios::binary); //abrimos el archivo de salida para escribir los resultados fusionados , va a ser el 
												 //archivo de salida con todos los resultados FUSIONADOS de cada rank
	for (int i = 0; i < num_ranks; ++i) //para cada rank 
	{
		std::ifstream infile(output_prefix_base + "_rank" + std::to_string(i) + extension,
				std::ios::in | std::ios::binary); // abrimos el archivo de salida de ese rank (sobre el que iteramos) 
		outfile << infile.rdbuf(); // copiamos el contenido del archivo de ese rank (sobre el q iteramos) al archivo de salida fusionado abierto antes
		infile.close(); //cerramos el archivo de entrada del rank actual (sobre el que iteramos)
	}
	outfile.close(); //por ultimo, una vez copiados todos los archivos de salida de cada rank, cerramos el archivo de salida FUSIONADO
}


int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv); // Inicializa MPI 

	int rank, size; //declaramos la variable rank para el id del proceso y size para el numero de procesos (rank representara el id de cada proceso (rank) el 0,1,...) 
			// size representara el numero total de procesos (ranks) que se estan ejecutando en paralelo num_ranks = size
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); //obtenemos el id del proceso actual (rank) y lo guardamos en la variable rank
	MPI_Comm_size(MPI_COMM_WORLD, &size); //obtenemos el numero total de procesos (ranks) que se estan ejecutando en paralelo y lo guardamos en la variable size

	auto config_opt = GanonClassify::CommandLineParser::parse(argc, argv); //parseamos los argumentos de la linea de comandos para obtener la configuracion de GanonClassify
	if (!config_opt.has_value()) { //si no se pudo parsear la configuracion, significa que hubo un error en los argumentos de la linea de comandos
		MPI_Finalize(); //finalizamos MPI antes de salir
		return argc == 1 ? EXIT_FAILURE : EXIT_SUCCESS; //esto es tal cual el codigo original solo que añado la finalizacion de MPI antes de salir
	}

	auto config = config_opt.value(); //obtenemos la configuracion de GanonClassify a partir del resultado del parseo de los argumentos de la linea de comandos

	std::string input_file = config.input_file; //obtenemos el archivo de entrada (el FASTQ o FASTA) de la config que se ha obtenido en el parseo de la linea de comandos de los 
						    //argumentos pasados al programa
	std::string output_prefix_base = config.output_prefix; //obtenemos el prejijo que se quiere usar para los archivos de salida obtenido en el parseo de la linea de comandos
	std::string output_prefix_rank = output_prefix_base + "_rank" + std::to_string(rank); //concatenamos el prefijo de la salida con el _rankn donde n es el id del rank
											      //cada rank un archivo de salida que luego el rank 0 se encargara de fusionar
	std::string chunked_input = "input_chunk_rank" + std::to_string(rank) + ".fq"; //nombre del archivo de entrada que va a usar cada ranks
	
	//OJO fq es abreviacion de FASTQ, como luego para procesar el archivo se usa seqan 3 que detecta el tipo de archivo por el contenido, no hay problemas con la extension

	// Solo el rank 0 divide el archivo de entrada
	if (rank == 0) {
		std::ifstream infile(input_file); //abrimos el archivo de entrada (el FASTQ ) para leerlo
		std::vector<std::ofstream> out_chunks(size); //vector de archivos de salida (chunks) que va a crear el rank 0, cada uno para cada rank
		for (int i = 0; i < size; ++i) //para cada rank
			out_chunks[i].open("input_chunk_rank" + std::to_string(i) + ".fq"); //abrimos el archivo de salida (chunk) para ese rank, cada rank tendra su propio archivo de entrada (chunk) que luego procesara

		std::string line; //variable para leer las lineas del archivo de entrada con todas las secuencias a clasificar
		int read_count = 0; //contador de lecturas (reads) que se han leido del archivo de entrada
		int line_idx = 0; //contador de lineas leidas del archivo de entrada, para saber cuando hemos leido 4 lineas (una secuencia completa en formato FASTQ)
		//Mas info del formato fastq en https://es.wikipedia.org/wiki/Formato_FASTQ
		//
		//mirar cuantas lineas tiene el archivo de entrada, y repartirlas entre los num_ranks
		//cada rank va a leer 4 lineas del archivo de entrada (una secuencia completa en formato FASTQ) y escribirla en su propio archivo de salida (chunk) 
		int lineasArchivoEntrada = std::count(std::istreambuf_iterator<char>(infile), std::istreambuf_iterator<char>(), '\n'); //contamos las lineas del archivo de entrada 
		infile.clear(); //limpiamos el estado del stream para poder volver a leerlo desde el principio
		infile.seekg(0); //volvemos al principio del archivo de entrada para poder leerlo desde el principio
		int lineasReales = lineasArchivoEntrada / 4; //calculamos el numero de secuencias completas en el archivo de entrada, ya que cada secuencia ocupa 4 lineas en formato FASTQ
		//toca calcular numero de lineas reales (de secuencias) que le toca a cada rank 
		if (lineasReales < size) { //si el numero de secuencias es menor que el numero de ranks, entonces cada rank va a leer una secuencia completa (4 lineas)
			size = lineasReales; //ajustamos el numero de ranks al numero de secuencias reales
		}
		//variable de lineas por rank 
		int lineasRealesPorRank = 0;
		
		if (lineasReales % size != 0) {
			//si el numero de lineas no es divisible por el numero de ranks, entonces un rank hara mas trabajo que los demas
			std::cerr << "Warning: The number of sequences is not evenly divisible by the number of ranks. Some ranks will process more sequences than others.\n";
			int aux = lineasReales % size; //guardamos el resto de la division para saber cuantas secuencias mas le toca a un ranks
			//el rank 0 procesara las que le toca mas el resto 
			if(rank == 0) {
				lineasRealesPorRank = lineasReales / size + aux; //el rank 0 procesara las secuencias que le tocan mas el resto
			} else {
				lineasRealesPorRank = lineasReales / size; //los demas ranks procesaran las secuencias que les tocan
			}
		}else {
			lineasRealesPorRank = lineasReales / size; //si el numero de lineas es divisible por el numero de ranks, entonces cada rank va a leer el mismo numero de secuencias
		}

		while (std::getline(infile, line)) { //leemos el archivo de entrada linea por linea
			out_chunks[read_count / lineasRealesPorRank ] << line << "\n"; //escribimos la linea en el archivo de salida (chunk) correspondiente al rank que le toca
			++line_idx; //incrementamos el contador de lineas leidas
			if (line_idx == 4) { //si hemos leido 4 lineas (una secuencia completa en formato FASTQ)
				line_idx = 0; //reiniciamos el contador de lineas leidas
				++read_count; //incrementamos el contador de lecturas (reads) leidas
			}
		}

		for (auto& f : out_chunks) // cerramos todos los archivos de salida (chunks) que hemos abierto para cada rank
			f.close();
	}

	MPI_Barrier(MPI_COMM_WORLD); // Esperamos a que todos los ranks terminen de dividir el archivo de entrada antes de continuar

	// Cada rank usa su chunk y su propio output_prefix
	config.input_file = chunked_input;
	config.output_prefix = output_prefix_rank;
	config.threads = 1;  // para no interferir con otros procesos MPI, de momento se va a hacer uso de 1 solo thread de la implementacion mutlhilo de ganon de hilos de C++
	config.quiet = true;

	bool result = GanonClassify::run(std::move(config)); // Ejecutamos GanonClassify con la configuracion obtenida y el chunk de entrada correspondiente a cada rank

	MPI_Barrier(MPI_COMM_WORLD); // Esperamos a que todos los ranks terminen de procesar su chunk antes de continuar

	if (rank == 0) // Solo el rank 0 fusiona los resultados
	{
		std::cout << "Fusionando resultados...\n"; //mensaje de inicio de fusionado de resultados
		merge_files(output_prefix_base, ".rep", size);
		merge_files(output_prefix_base, ".lca", size); //comprobar si hay arhivos lca ??????
		merge_files(output_prefix_base, ".all", size);
		merge_files(output_prefix_base, ".tre", size);

		//merge_files(output_prefix_base, ".unc", size);
		std::cout << "Listo. Todos los archivos parciales de los ranks fueron fusionados\n";
	}

	MPI_Finalize(); // Finaliza MPI antes de salir
	return result ? EXIT_SUCCESS : EXIT_FAILURE;
}
// Este codigo es una primera aproximacion a la paralelizacion de GanonClassify usando MPI

#=
  Script: funcionesLecturaEscritura.jl
  Autor: Adrián Martín
  Última revisión: 15/03/2021
  Coms: este script implementa funciones que abstraen la lectura 
        de ficheros en formato FASTA y CSV, así como la escritura 
        de estructuras de datos matriciales en ellos
=#

using DelimitedFiles    # parser de ficheros estructurados en filas y columnas
using BioSequences      # incluye parser de formato FASTA

#**********************************************************************************************************************
# TIPO DE DATO PARA FACILITAR EL MANEJO REQUERIDO DE LAS SECUENCIAS
mutable struct Secuencia
        identificador::String
        secuencia::DNASequence
end

#**********************************************************************************************************************
# FUNCIONES PARA LECTURA DE FICHEROS

# Lee el fichero CSV indicado en <<rutaFicheroCSV>> y devuelve la matriz resultante de parsearlo como Float
function leerCSV(rutaFicheroCSV)
        try
	        global matriz = readdlm(rutaFicheroCSV, ',', Float64, '\n')
        catch
	        println("ERROR: no se ha podido leer el fichero ", rutaFicheroCSV)
	        exit(1) # terminar la ejecución de Julia
        end

        return matriz
end

# Lee el fichero FASTA indicado en <<rutaFicheroFASTA>> y devuelve la matriz que tiene por filas las secuencias almacenadas en él
function leerFASTA(rutaFicheroFASTA)
        try
	        # De esta forma, el stream del fichero solo está abierto hasta finalizar el ámbito del bloque do..end
                open(FASTA.Reader, rutaFicheroFASTA) do reader
                        record = read(reader) # lectura de la primera secuencia
                        global secuencias = [Secuencia(FASTA.identifier(record), FASTA.sequence(record))] # al definirlo como global, es accesible desde el ámbito raiz
                        for record in reader  # lectura del resto de secuencias
                                secuencias = [secuencias ; Secuencia(FASTA.identifier(record), FASTA.sequence(record))] # añade una fila en cada iteración
                        end
                end
        catch
	        println("ERROR: no se ha podido leer el fichero ", rutaFicheroFASTA)
	        exit(1) # terminar la ejecución de Julia
        end

        return secuencias
end

#**********************************************************************************************************************
# FUNCIONES PARA ESCRITURA EN FICHEROS

# Escribe la estructura de datos <<matriz>> en formato CSV en un fichero de nombre <<rutaFicheroCSV>> 
function escribirCSV(rutaFicheroCSV, matriz)
        try
                writedlm(rutaFicheroCSV, matriz, ',')
        catch
                println("ERROR: no se han podido escribir el fichero ", rutaFicheroCSV)
                exit(1) # terminar la ejecución de Julia
        end

        return nothing
end

# Escribe las secuencias dadas en <<secuencias>> en formato FASTA en un fichero de nombre <<rutaFicheroFASTA>>
function escribirFASTA(rutaFicheroFASTA, secuencias)
        try
	        # De esta forma, el stream del fichero solo está abierto hasta finalizar el ámbito del bloque do..end
                open(FASTA.Writer, rutaFicheroFASTA) do writer
                        for secuencia in secuencias
                                record = FASTA.Record(secuencia.identificador, secuencia.secuencia)
                                write(writer, record)
                        end
                end
        catch
	        println("ERROR: no se ha podido escribir el fichero ", rutaFicheroFASTA)
	        exit(1) # terminar la ejecución de Julia
        end

        return nothing
end

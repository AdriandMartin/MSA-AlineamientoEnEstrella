#!/usr/bin/julia -O3
#=
  Script: calcularSP.jl
  Autor: Adrián Martín
  Última revisión: 16/03/2021
  Coms: este script calcula la puntuación por parejas, SP, de un 
        alineamiento a partir de su matriz de puntuación
        La interpretación del resultado dependerá de la métrica 
        usada para calcular los elementos de la matriz dada  
=#

include("lib/funcionesLecturaEscritura.jl")

# Comprobación de argumentos de entrada (debe aportarse uno, que es la ruta al fichero con la matriz de puntuaciones)
if size(ARGS,1) != 1
	println("ERROR: usar como \"", PROGRAM_FILE, " ruta/al/ficheroConMatrizPuntuaciones\"")
	exit(1)
end 

# Lectura del fichero y asignación de su contenido a una matriz de float
matriz = leerCSV(ARGS[1])

M = size(matriz,1) # número de secuencias representadas en la matriz

# Obtención de la puntuación por parejas
puntuacion = sum(matriz) / (M * (M-1)) # equivale a "(sum(matriz) / 2) / ((M * (M-1)) / 2)", es decir, a la normalización de los elementos de uno de los triángulos simétricos de la matriz

println("SP: ", puntuacion)

exit(0)


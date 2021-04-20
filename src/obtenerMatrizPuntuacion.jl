#!/usr/bin/julia -O3
#=
  Script: obtenerMatrizPuntuacion.jl
  Autor: Adrián Martín
  Última revisión: 15/03/2021
  Coms: este script calcula la puntuación por parejas de un 
        alineamiento a partir de sus secuencias, generando como 
        resultado la matriz de puntuación correspondiente
    Se consideran los valores 'N' como errores de secuenciación, 
        y se interpretan como 'cualquier base'
    El sistema de puntuación que utiliza implementa penalización 
        afín, pero la puntuación para cada tipo de columna es 
        configurable
  Nota: para activar el log de debug, ejecutar en shell antes 
        de lanzar este script lo siguiente (añade la variable 
        JULIA_DEBUG a las variables de entorno con ese valor)
        export JULIA_DEBUG="all"
        para dejar de ver el log de debug, darle valor ""
=#

using BioSequences      # definición de dna_str
using ArgParse          # parser de argumentos por línea de comandos
using Combinatorics     # incluye funciones para generar todas las combinaciones y permutaciones de un conjunto

include("lib/funcionesLecturaEscritura.jl")

#**********************************************************************************************************************
# TIPO DE DATO PARA FACILITAR LA REPRESENTACIÓN DEL SISTEMA DE PUNTUACIÓN
# Aviso: se usa este struct y no BioAlignments.AffineGapScoreModel dada a la complejidad innecesaria de este último para lo que se va a hacer
struct Puntuacion
    match::Int
    mismatch::Int
    gap_open::Int
    gap_extend::Int
end

#=**********************************************************************************************************************=#
# FUNCIONES

# Función que define y realiza el parsing de los argumentos por línea de comandos
function parse_commandline()
    s = ArgParseSettings()

    s.description = "Este script calcula la puntuación por parejas de un alineamiento a partir del mismo, generando como resultado una matriz de puntuación (simétrica y con diagonal cero)"
    s.epilog = "El sistema de puntuación que utiliza implementa penalización afín, pero la puntuación para cada tipo de columna es configurable"

    @add_arg_table! s begin
        "--match"
            help = "puntuación en caso de coincidencia"
            arg_type = Int
            default = 1
        "--mismatch"
            help = "puntuación en caso de no coincidencia"
            arg_type = Int
            default = -3
        "--gap_open"
            help = "puntuación de cada primer gap"
            arg_type = Int
            default = -5
        "--gap_extend"
            help = "puntuación de los sucesivos gaps"
            arg_type = Int
            default = -2
        "INPUT_FILE"
            help = "fichero FASTA con el alineamiento a evaluar"
            arg_type = String
            action = :store_arg # guardar como variable su valor
            required = true
        "-o","--output"
            help = "fichero en el que volcar la matriz de puntuación entre pares"
            arg_type = String
            dest_name = "OUTPUT_FILE"
            default = "score_matrix"
    end

    return parse_args(s)
end

# Enumeración con los tres posibles estados durante la comparación de las bases de las secuencias
@enum Estado GAP_DETECTADO_SEQA GAP_DETECTADO_SEQB ESPERANDO_GAP

# Evalúa las secuencias seqA y seqB según el sistema de puntuación dado
# Pre: se asume que seqA y seqB están alineadas
function puntuar(seqA::DNASequence, seqB::DNASequence, puntuacion::Puntuacion)::Int
    estado = ESPERANDO_GAP
    score = 0
    for i in 1:length(seqA) # length(seqA) = length(seqB)
        if seqA[i] == seqB[i]
            if seqA[i] != DNA_Gap && seqA[i] != DNA_N 
                score += puntuacion.match
                estado = ESPERANDO_GAP
            # else score += 0 # caso especial, posiciones no se cuentan
            end
        else # seqA[i] != seqB[i]
            if (seqA[i] == DNA_Gap) # seqB[i] == {A,C,G,T,N}
                if estado == GAP_DETECTADO_SEQA
                    score += puntuacion.gap_extend
                else
                    score += puntuacion.gap_open
                end
                estado = GAP_DETECTADO_SEQA
            elseif (seqB[i] == DNA_Gap) # seqA[i] == {A,C,G,T,N}
                if estado == GAP_DETECTADO_SEQB
                    score += puntuacion.gap_extend
                else
                    score += puntuacion.gap_open
                end
                estado = GAP_DETECTADO_SEQB
            elseif (seqA[i] == DNA_N) || (seqB[i] == DNA_N) # no hay gaps, pero una de las bases es N (no las dos, porque son distintas entre sí)
                score += puntuacion.match
                estado = ESPERANDO_GAP
            else # no hay ni gaps ni Ns, son letras A,C,G o T distintas entre sí
                score += puntuacion.mismatch
                estado = ESPERANDO_GAP
            end 
        end
        @debug println("(",seqA[i],",",seqB[i],") -> ",score)
    end
    return score
end

#=**********************************************************************************************************************=#
# MAIN

# Parsing de argumentos de línea de comandos
parsed_args = parse_commandline()
alignment_file = parsed_args["INPUT_FILE"]
output_file = parsed_args["OUTPUT_FILE"]
match_score = parsed_args["match"]
mismatch_score = parsed_args["mismatch"]
gap_open_score = parsed_args["gap_open"]
gap_extend_score = parsed_args["gap_extend"]

# Crear sistema de puntuación
puntuacion = Puntuacion(match_score, mismatch_score, gap_open_score, gap_extend_score)

# Leer secuencias
secuencias = leerFASTA(alignment_file)
M = size(secuencias,1)

# Comprobación de que todas las secuencias tienen la misma longitud (están alineadas)
longitud = length(secuencias[1].secuencia)
for i in 2:M
    @assert longitud == length(secuencias[i].secuencia)
end

# Construcción de matriz de puntuación
matriz = zeros(M,M)
Threads.@threads for (i,j) = collect(combinations(1:M,2))
    matriz[i,j] = matriz[j,i] = puntuar(secuencias[i].secuencia, secuencias[j].secuencia, puntuacion)
end

# Volcar matriz de score en el fichero indicado
escribirCSV(output_file, matriz)

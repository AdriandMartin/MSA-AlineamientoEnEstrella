#!/usr/bin/julia -O3
#=
  Script: detectarB117.jl
  Autor: Adrián Martín
  Última revisión: 16/03/2021
  Coms: este script busca mutaciones correspondientes a la variante B.1.1.7 o 
        británica sobre secuencias de SARS-COV-2 alineadas
    Nótese que las secuencias dadas como entrada deben estar alineadas, y la 
        primera de ellas debe ser la secuencia de referencia de SARS-COV-2
  Posible mejora: generar una tabla con los resultados en lugar de solo 
                  mostrarlos por pantalla
=#

using BioSequences      # definición de dna_str
using ArgParse          # parser de argumentos por línea de comandos

include("lib/funcionesLecturaEscritura.jl")

#**********************************************************************************************************************
# FUNCIONES

# Función que define y realiza el parsing de los argumentos por línea de comandos
function parse_commandline()
  s = ArgParseSettings()

  s.description = "Este script busca mutaciones correspondientes a la variante B.1.1.7 o británica sobre secuencias de SARS-COV-2 alineadas"
  s.epilog = "Nótese que las secuencias dadas como entrada deben estar alineadas, y la primera de ellas debe ser la secuencia de referencia de SARS-COV-2"

  @add_arg_table! s begin
      "--ORF1ab"
          help = "buscar mutaciones en el gen ORF1ab"
          action = :store_true
      "--Spike"
          help = "buscar mutaciones en el gen Spike"
          action = :store_true
      "--ORF8"
          help = "buscar mutaciones en el gen ORF8"
          action = :store_true
      "--N"
          help = "buscar mutaciones en el gen N"
          action = :store_true
      "-a", "--all"
          help = "buscar mutaciones en todos los genes: ORF1ab, Spike, ORF8 y N"
          action = :store_true
      "INPUT_FILE"
          help = "fichero FASTA con las secuencias alineadas"
          arg_type = String
          action = :store_arg # guardar como variable su valor
          required = true
  end

  return parse_args(s)
end

# Devuelve el mapeo de posiciones en la secuencia dada con respecto a su equivalente sin gaps
function mapearPosiciones(seq::DNASequence)
    rseq2seq = Dict{Int,Int}() # a cada posición en la secuencia sin gaps se le asocia la posición en la secuencia con gaps
    i = 1    # índice sobre la secuencia sin gaps
    real = 1 # índice sobre la secuencia dada (con gaps)
    while real <= length(seq)
        if seq[real] != DNA_Gap # seq[real] es una base y aparece en la secuencia sin gaps
            rseq2seq[i] = real # se asocia a esa posición sin gaps la real en la secuencia
            real += 1
            i += 1
        else # seq[real] es un gap, luego esa posición no tiene asociada en la secuencia sin gaps
            real += 1
        end
    end
    return rseq2seq
end

# Buscar mutación del gen ORF1ab en las secuencias dadas
function buscarORF1ab(secuencias::Array{Secuencia}, ref2seq::Dict{Int,Int})
    println("Buscando mutaciones en el gen ORF1ab...")
    # Buscar C3267T
    println("> Sustitución C → T en posición 3267")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[3267]] == DNA_T)
    end
    # Buscar C5388A
    println("> Sustitución C → A en posición 5388")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[5388]] == DNA_A)
    end
    # Buscar T6954C
    println("> Sustitución T → C en posición 6954")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[6954]] == DNA_C)
    end
    # Buscar borrado (gaps) de 11288 a 11296
    println("> Borrado de posición 11288 a la 11296")
    for s in secuencias
        borrado = [s.secuencia[ref2seq[i]] == DNA_Gap for i in 11288:11296]
        println("\t", s.identificador, " -> ", all(borrado))
    end
end

# Buscar mutación del gen Spike en las secuencias dadas
function buscarSpike(secuencias::Array{Secuencia}, ref2seq::Dict{Int,Int})
    println("Buscando mutaciones en el gen Spike...")
    # Buscar borrado (gaps) de 21765 a 21770
    println("> Borrado de posición 21765 a la 21770")
    for s in secuencias
        borrado = [s.secuencia[ref2seq[i]] == DNA_Gap for i in 21765:21770]
        println("\t", s.identificador, " -> ", all(borrado))
    end
    # Buscar borrado (gaps) de 21991 a 21993
    println("> Borrado de posición 21991 a la 21993")
    for s in secuencias
        borrado = [s.secuencia[ref2seq[i]] == DNA_Gap for i in 21991:21993]
        println("\t", s.identificador, " -> ", all(borrado))
    end
    # Buscar A23063T
    println("> Sustitución A → T en posición 23063")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[23063]] == DNA_T)
    end
    # Buscar C23271A
    println("> Sustitución C → A en posición 23271")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[23271]] == DNA_A)
    end
    # Buscar C23604A
    println("> Sustitución C → A en posición 23604")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[23604]] == DNA_A)
    end
    # Buscar C23709T
    println("> Sustitución C → T en posición 23709")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[23709]] == DNA_T)
    end
    # Buscar T24506G
    println("> Sustitución T → G en posición 24506")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[24506]] == DNA_G)
    end
    # Buscar G24914C
    println("> Sustitución G → C en posición 24914")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[24914]] == DNA_C)
    end
end

# Buscar mutación del gen ORF8 en las secuencias dadas
function buscarORF8(secuencias::Array{Secuencia}, ref2seq::Dict{Int,Int})
    println("Buscando mutaciones en el gen ORF8...")
    # Buscar C27972T
    println("> Sustitución C → T en posición 27972")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[27972]] == DNA_T)
    end
    # Buscar G28048T
    println("> Sustitución G → T en posición 28048")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[28048]] == DNA_T)
    end
    # Buscar A28111G
    println("> Sustitución A → G en posición 28111")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[28111]] == DNA_G)
    end
end

# Buscar mutación del gen N en las secuencias dadas
function buscarN(secuencias::Array{Secuencia}, ref2seq::Dict{Int,Int})
    println("Buscando mutaciones en el gen N...")
    # Buscar G28280C
    println("> Sustitución G → C en posición 28280")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[28280]] == DNA_C)
    end
    # Buscar A28281T
    println("> Sustitución A → T en posición 28281")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[28281]] == DNA_T)
    end
    # Buscar T28282A
    println("> Sustitución T → A en posición 28282")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[28282]] == DNA_A)
    end
    # Buscar C28977T
    println("> Sustitución C → T en posición 28977")
    for s in secuencias
        println("\t", s.identificador, " -> ", s.secuencia[ref2seq[28977]] == DNA_T)
    end
end

#=**********************************************************************************************************************=#
# MAIN

# Parsing de argumentos de línea de comandos
parsed_args = parse_commandline()
input_file = parsed_args["INPUT_FILE"]

# Leer secuencias
secuencias = leerFASTA(input_file)
M = size(secuencias,1)

# Comprobación de que todas las secuencias tienen la misma longitud (están alineadas)
longitud = length(secuencias[1].secuencia)
for i in 2:M
    @assert longitud == length(secuencias[i].secuencia)
end

# Extracción de la cadena de DNA correspondiente a la secuencia de referencia
refSeq = secuencias[1].secuencia

# Mapear posiciones de la secuencia de referencia original en la secuencia con gaps alineada
ref2seq = mapearPosiciones(refSeq)

# Buscar las mutaciones en los genes indicados por línea de comandos
if parsed_args["all"] || parsed_args["ORF1ab"]
    buscarORF1ab(secuencias[2:end], ref2seq)
end

if parsed_args["all"] || parsed_args["Spike"]
    buscarSpike(secuencias[2:end], ref2seq)
end

if parsed_args["all"] || parsed_args["ORF8"]
    buscarORF8(secuencias[2:end], ref2seq)
end

if parsed_args["all"] || parsed_args["N"]
    buscarN(secuencias[2:end], ref2seq)
end

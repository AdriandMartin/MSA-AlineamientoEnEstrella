#!/usr/bin/julia -O3
#=
  Script: estrellaMSA.jl
  Autor: Adrián Martín
  Última revisión: 16/03/2021
  Coms: este script realiza el multialineamiento de las secuencias dadas, 
        usando para ello el método de la estrella, con centro otra secuencia 
        aportada
    El sistema de puntuación que utiliza implementa penalización afín, pero 
        la puntuación para cada tipo de columna es configurable
  Posible mejora: en caso de no aportarse una secuencia de referencia, 
                  tomar como centro la secuencia a menor distancia del 
                  resto entre las dadas
=#

using BioAlignments     # definición de pairalign y los TADs usados por éste
using BioSequences      # definición de dna_str
using ArgParse          # parser de argumentos por línea de comandos

include("lib/funcionesLecturaEscritura.jl")

#**********************************************************************************************************************
# TIPO DE DATO PARA FACILITAR EL MANEJO REQUERIDO DE LOS ALINEAMIENTOS
mutable struct Alineamiento
  referencia::Secuencia
  resto::Array{Secuencia}
  bp::Int # número de pares de bases

  añadir_gap::Function
  añadir_secuencias::Function

  #------------------------------------------
  # Función para añadir un gap en la misma posición de todas las secuencias del alineamiento
  # Pre: posicion es un valor entero entre 0 y la longitud de las secuencias del alineamiento
  # Nota: add_gap(aln,i) supone añadir un gap que ocupará la posición i+1 en las secuencias 
  #       del alineamiento aln una vez insertado
  function add_gap(self::Alineamiento, posicion::Int)
    @assert 0 <= posicion && posicion <= self.bp
    self.referencia.secuencia = self.referencia.secuencia[1:posicion] * dna"-" * self.referencia.secuencia[posicion+1:self.bp]
    for i = 1:size(self.resto,1)
      self.resto[i].secuencia = self.resto[i].secuencia[1:posicion] * dna"-" * self.resto[i].secuencia[posicion+1:self.bp]
    end
    self.bp += 1
  end

  # Función para añadir un vector <<sequences>> de secuencias al vector <<self.resto>> de secuencias
  # del alineamiento
  # Pre: las secuencias a añadir al alineamiento deben estar alineadas con las que éste ya almacena
  function add_sequences(self::Alineamiento, sequences::Array{Secuencia})
    for sequence in sequences
      @assert length(sequence.secuencia) == self.bp
    end
    self.resto = [self.resto; sequences]
  end

  #------------------------------------------
  # Constructor
  function Alineamiento(refSeq::Secuencia, seq::Secuencia)
    longitud = length(refSeq.secuencia)
    @assert longitud == length(seq.secuencia)
    new(refSeq, [seq], longitud, add_gap, add_sequences)
  end

end

#**********************************************************************************************************************
# FUNCIONES

# Función que define y realiza el parsing de los argumentos por línea de comandos
function parse_commandline()
  s = ArgParseSettings()

  s.description = "Este script realiza el multialineamiento de las secuencias dadas, usando para ello el método de la estrella, con centro otra secuencia aportada"
  s.epilog = "NOTA: Julia por defecto ejecuta todo en un único hilo, pero este script está implementado de forma que permite paralelizar ciertas instrucciones iterativas. Para indicar el número de hilos a utilizar, ejecutar antes en el shell 'export JULIA_NUM_THREADS=N', siendo N un valor natural menor o igual que el número de núcleos físicos de la máquina (por experiencia, se recomienda que sea menor)"

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
      "-r","--ref-seq"
          help = "fichero FASTA con la secuencia de referencia a utilizar como centro de la estrella"
          arg_type = String
          dest_name = "REFSEQ_FILE"
          action = :store_arg # guardar como variable su valor
          required = true
      "INPUT_FILE"
          help = "fichero FASTA con las secuencias a alinear"
          arg_type = String
          action = :store_arg # guardar como variable su valor
          required = true
      "-o","--output"
          help = "fichero FASTA en el que volcar el alineamiento"
          arg_type = String
          dest_name = "OUTPUT_FILE"
          default = "alignment.fasta"
  end

  return parse_args(s)
end

# Devuelve una instancia de Alineamiento con el resultado de alinear dos secuencias dadas usando el sistema de puntuación indicado
function alinear(refSeq::Secuencia, seq::Secuencia, puntuacion::AffineGapScoreModel)::Alineamiento
  resultado = pairalign(GlobalAlignment(), refSeq.secuencia, seq.secuencia, puntuacion)
  aln = alignment(resultado)
  refSeqAln = DNASequence([x for (x, _) in aln])
  seqAln = DNASequence([x for (_, x) in aln])
  alineamiento = Alineamiento(Secuencia(refSeq.identificador,
                                        refSeqAln), 
                              Secuencia(seq.identificador,
                                        seqAln)
  )
  return alineamiento
end

# Devuelve el alineamiento resultante de combinar los dos dados usando "Once a gap - always a gap"
function combinar(alnI::Alineamiento, alnJ::Alineamiento)::Alineamiento
  i = 1 # índice con el que recorrer alnI y alnJ
  # Recorrer ambos alineamientos añadiendo los gaps necesarios
  while (i <= alnI.bp) && (i <= alnJ.bp) # i <= min(alnI.bp,alnJ.bp)
    # Comparación de las bases apuntadas por el cursor
    if alnI.referencia.secuencia[i] != alnJ.referencia.secuencia[i]
      if alnI.referencia.secuencia[i] == DNA_Gap
        # Se necesita añadir un gap en alnJ
        alnJ.añadir_gap(alnJ, i-1)
      elseif alnJ.referencia.secuencia[i] == DNA_Gap # no queda otra, porque son la misma secuencia
        # Se necesita añadir un gap en alnI
        alnI.añadir_gap(alnI, i-1)
      end
    end
    # Se avanza el cursor
    i += 1
  end
  # Arreglar situación en la que una secuencia tenía gaps al final
  if alnI.bp != alnJ.bp
    if alnI.bp < alnJ.bp
      while alnI.bp != alnJ.bp
        alnI.añadir_gap(alnI, i-1)
        i += 1
      end
    else # alnI.bp > alnJ.bp
      while alnI.bp != alnJ.bp
        alnJ.añadir_gap(alnJ, i-1)
        i += 1
      end
    end
  end
  # Comprobar alineamiento y unirlos en una sola instancia
  @assert alnI.referencia.secuencia == alnJ.referencia.secuencia
  alnI.añadir_secuencias(alnI,alnJ.resto)
  return alnI
end

# Devuelve la lista de secuencias dadas alineadas aplicando el método de la estrella, usando como centro la secuencia <<referencia>>
function estrella(referencia::Secuencia, secuencias::Array{Secuencia})::Array{Secuencia}
  # Obtención de los alineamientos por pares
  alineamientos = Array{Alineamiento,1}(undef,size(secuencias,1))
  Threads.@threads for i in 1:size(secuencias,1) # de esta forma, se reparten las iteraciones de forma equitativa entre los JULIA_NUM_THREADS threads configurados
    alineamientos[i] = alinear(referencia,secuencias[i],puntuacion)
  end

  # Unión de los alineamientos por pares en uno usando "Once a gap - always a gap"
  alineamiento = alineamientos[1]
  for i in 2:size(alineamientos,1)
    alineamiento = combinar(alineamiento, alineamientos[i])
  end

  return [[alineamiento.referencia] ; alineamiento.resto]
end

#=**********************************************************************************************************************=#
# MAIN

# Parsing de argumentos de línea de comandos
parsed_args = parse_commandline()
input_file = parsed_args["INPUT_FILE"]
output_file = parsed_args["OUTPUT_FILE"]
refseq_file = parsed_args["REFSEQ_FILE"]
match_score = parsed_args["match"]
mismatch_score = parsed_args["mismatch"]
gap_open_score = parsed_args["gap_open"]
gap_extend_score = parsed_args["gap_extend"]

# Métrica de puntuación (maximizar) 
puntuacion = AffineGapScoreModel(
    match=match_score,
    mismatch=mismatch_score,
    gap_open=gap_open_score,
    gap_extend=gap_extend_score
)

# Leer secuencias
referencia = leerFASTA(refseq_file)[1] # quedarse con la secuencia como tal
secuencias = leerFASTA(input_file)

# Obtener las secuencias alineadas usando el método de la estrella
println("Alineando secuencias de '",input_file,"' con '",refseq_file,"' como secuencia de referencia, usando ",Threads.nthreads()," hilos...")
texe = @elapsed secuenciasAlineadas = estrella(referencia, secuencias)
println("Multialineamiento completado en ",texe," segundos")

# Escribir multi-alineamiento resultante en el fichero indicado
println("Escribiendo resultado en fichero '",output_file,"'")
escribirFASTA(output_file, secuenciasAlineadas)

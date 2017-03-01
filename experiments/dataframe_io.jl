#using DataFrames
export read_header_lines, write_data_frame, write_dataframe, read_headers

function read_header_lines( stream::IO )
  result = String[]
  line = readline( stream )
  while line[1] == '#'
    push!( result, line )
    line = readline( stream )
  end
  result
end

function add_quotes( string_array::Array{String,1} )
  map(n->"\"$n\"",string_array)
end

function write_data_frame( stream::IO, header_lines::Vector{String}, df::DataFrame )
  for h in header_lines
    if h[1] == '#' && h[end] == '\n'
      write(stream, h )
    else
      write(stream, '#', h, "\n" )
    end
  end 
  write(stream, join( add_quotes(map(string,names(df))), "," ), "\n" )
  for i = 1:size(df)[1]
    write(stream, join( Any[ df[j][i] for j = 1:size(df)[2]], "," ), "\n" )
  end
  close(stream)
end
  
function write_dataframe( filename::String, header_lines::Vector{String}, df::DataFrame )
  str = open(filename,"w")
  write_data_frame( str, header_lines, df )
  close(str)
end

function read_headers( filename::String)
  str = open(filename,"r")
  header_lines = read_header_lines( str )
  close(str)
  header_lines
end
  

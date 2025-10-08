module ConfigFile
"""
configfile.py - Human-readable text configuration file library 
Copyright 2010  Luke Campagnola
Distributed under MIT/X11 license. See license.txt for more infomation.

Used for reading and writing dictionary objects to a python-like configuration
file format. Data structures may be nested and contain any data type as long
as it can be converted to/from a string using repr and eval.
"""

using Formatting
using PyCall
pgc = pyimport("python/configfile")  # use python configuration reader... 

export test, pgtest, testparse1

# import re, os, sys, datetime
# import numpy
# from .pgcollections import OrderedDict
# from . import units
# from .python2_3 import asUnicode, basestring
# from .Qt import QtCore
# from .Point import Point
# from .colormap import ColorMap
# GLOBAL_PATH = None # so not thread safe.



# function writeConfigFile(data, fname):
#     s = genString(data)
#     fd = open(fname, "w")
#     fd.write(s)
#     fd.close()
#
    function readConfigFile(fname)
        #os.chdir(newDir)  ## bad.
        fd = open(fname)
        s = read(fd, String)
        close(fd)
        s = replace(s, "\r\n" => "\n")
        s = replace(s, "\r" => "\n")
        data = parseString(s)
        return data
    end
    # function appendConfigFile(data, fname)
    #     s = genString(data)
    #     open(fname, "a") do fd
    #         fd.write(s)
    #     end
    # end

    function genString(data, indent="")
        s = ""
        for k in data
            sk = str(k)
            if length(sk) == 0
                println("genstring: data =", data)
                println("Blank dict keys are not allowed (see data above).")
                return nothing
            end
            if (sk[1] == " ") | ":" in sk
                println("genstring: data = ", data)

                println("Dict keys must not contain ':' or start with spaces [offending key is ",  sk, "]")
                return nothing
            end
            if isinstance(data[k], dict)
                s = s * indent * sk * ":\n"
                s = s * genString(data[k], indent * "    ")
            else
                s = indent * sk * ": " * repr(data[k]) * "\n"
            end
        end
        return s
    end

    function printlines(lines)
        n = 1
        # println(length(lines))
        for line in lines
            println(n, ": ", line) # "{0:%3d}: {1:%64s}\n" % n, line)
            n += 1
        end
    end
    
    function removecomments(lines)
        out_lines = "" 
        # clean up the lines so that they only have non-commented lines
        # and remove the comment text
        # println("\nRML")
        for line in lines
            cleanline = strip(line, [' ', '\t'])
            # println("   rml: # length of cleaned line: ", length(cleanline))
            if length(cleanline) == 0
                continue
                end
            # println("   rml:: cleanline: ", cleanline)
            if startswith(cleanline,'#')
                continue
            end
            if occursin(r"\s*#", line)
                # println("   rml: line has a following comment")
                f = match(r"#\S*", line)
                line = line[1:f.offset-1]
                # println("   rml: trimmed line: ", line)
            end
            if occursin(r"\S", line) | occursin(r"\s*#", line)
                line = rstrip(line) # clean trailing spaces
                # println("   rml: accept line with: <", line, ">\n")
                #push!(outlines, line)  # build to output
                out_lines *= line * "\n"
            end
        end
        return split(out_lines, "\n")
    end
    
    function measureIndent(s)
        n = 0
        while (n < length(s)) & occursin(s[n+1], " \t")
            n += 1
        end
        return n
    end  
    
    function parseTuple(v)
        println("V: ", v)
        println(" v1 : ", v[1], "  vend: ", v[end])
        if occursin(v[1], "[") & occursin(v[end], "]")  # recurse with nested array/tuple
            v = v[2:end-1]  # remove delimiters
            v = [parseTuple(v)]
        end
        if occursin(v[1], "(") & occursin(v[end], ")")  # recurse with nested array/tuple
            v = v[2:end-1]  # remove delimiters
            v = (parseTuple(v))
        end
        if ! isa(v, String)
            return v
        end
        println("parsing tuple of: ", v, " type: ", typeof(v))
        if occursin(",", v)  # commas separating multiple elements
            v=replace(v, " " => "")  # remove spaces
        else
            v=replace(v, " " => ",") # no commas, so insert instead of spaces...
        end
        println("v = ", v)
        # look for single-quoted strings first
        f = findall(r"\'\S+\'", v) #| occursin("\"" in v)
        val = []
        println("f: ", f)
        if f != nothing
            println("found something", f)
            
            for x in f
                println(x)
                ns = replace(v[f[x]], "'" => "\"")
                print("ns: ", ns)
                push!(val, ns)
            end
            println(val)
            println(maximum(f[end])+1)
            v = v[maximum(f[end])+1:end]
        end

        if occursin(".", v) | occursin("e", v) | occursin("E", v)
            vd = parse.(Float64, split(v, ","))
            push!(val, vd)
        else
            vd = parse.(Int64, split(v, ","))
            push!(val, vd)
        end
        return val
    end    
    
    function parseString(lines_in::String; start=1)
        """
        start is the line we start at - allows recursion into nested keys
        
        """
        data = Dict{Any, Any}()
        lines = split(lines_in, "\n")  # then split along the newlines
        lines = removecomments(lines)
        # println("Cleaned: ")
        # printlines(lines)
        indent = measureIndent(lines[start])  # get curreent indentation
        ln = start
        line = ""  # keep line in scope
        # try
            while true
                # println(ln, " ln ", typeof(ln))
                # println(length(lines), " lines ", typeof(length(lines)))
                if ln >= length(lines)
                    # println("Reached end of file at line: ", ln)
                    break
                end
                # println("continuing")
                line = lines[ln]
                if length(line) == 0  # skip blank lines
                    ln += 1
                    continue
                end
                # println("line: ", ln, ": ", line)

                ## Measure line indentation, make sure it is correct for this level
                lineInd = measureIndent(line)
                if lineInd < indent
                    # ln -= 1
                    break
                end
                # println("indent: ", lineInd, ", ", indent)
                if lineInd > indent
                    #print lineInd, indent
                    println("Indentation is incorrect. Expected ", indent, " got ", lineInd, " gn line ",  ln+1, " ",  line)
                    return ln, nothing
                end

                if !occursin(":", line)
                    println("Missing colon", ln+1, " ", line)
                    return ln, nothing
                end
                # println("parsing line to split")
                sp = split(line, ":")
                if length(sp) < 2
                    println("Missing name preceding colon", ln+1, line)
                    return nothing
                end
                # println("split ok, got: ", sp)
                k = strip(sp[1])
                v = strip(sp[2])
                # println("split into: ", k, " and ", v)
                if length(v) == 0  # key with no entry means new dict: try indented parsing
                    # println("   GOING DEEPER >>>>> at line: ", ln)
                    ln, val = parseString(lines_in, start=ln+1)
                    data[k] = val
                    continue
                end
                if (k[1] == "(") & (k[end] == ")")  ## If trimmedhe key looks like a tuple, try evaluating it.
                    # println("trying to parse key tuple")
                    try
                        k1 = parse.(Tuple, k) # eval(k, local)
                        if istype(k1, Tuple)
                            k = k1
                        end
                    catch
                        println("Failed to read tuple format of Key: ", k)
                        return ln, nothing
                    end
                end
                # println("    Evaluate the value")
                v = replace(v, "L" => "")  # remove "long" type markers on integers
                if occursin(v[1], "[(") & occursin(v[end], "])")  # is a tuple or an array? 
                    try
                        val = parseTuple(v)
                    # println("(nc) setting key k with value", k, " ", val)
                        data[k] = val
                    catch
                        val = Meta.parse(v) |> eval
                        data[k] = val
                    finally
                        println("failed to parse array of some sort: ", v)
                        return nothing
                    end
                    # println()
                elseif occursin(v[1], "\"") & occursin(v[end], "\"")  # quoted string
                    val = v[2:end-1]  # just get the characters
                    data[k] = val
                elseif occursin("'", v) #| occursin("\"" in v)
                    val = replace(v, "'" => "\"")[2:end] # v[2:end-1]
                    data[k] = val
                elseif occursin("False", v) | occursin("True", v)
                    v = lowercase(v)
                    val = parse(Bool, v)
                    data[k] = val
                elseif occursin(".", v) | occursin("e", v) | occursin("E", v) # float number
                    val = parse(Float64, v)
                    data[k] = val
                elseif (! occursin(",", v))  # try integer
                    val = parse(Int64, v)
                    data[k] = val
                elseif (ln >= length(lines)) | (measureIndent(lines[ln]) <= indent)
                        break # val = Dict()
                # elseif measureIndent(lines(ln) > indent)
                #         println( "Going deeper..", ln+1)
                #         (ln, val) = parseString(lines_in, start=ln)
                end
                 # println("(c) setting key", k, " with value", val)
                 data[k] => val
                # println("ln: ", ln, "  Data: ", data)
                ln += 1
            end # of while
            #print k, repr(val)
        # catch
        #     println("ParseError on  line: ", (ln)," ", line)
        # end
        # println("Returning  data = ", data)
        return ln, data
    end


    function testparse1()
        v = "(9, 3, 8)"
        val = parseTuple(v)
        println(v," of type: ", typeof(v)," becomes: ", val, " of type: ", typeof(val))
        
        v = "[1 3]"
        val = parseTuple(v)
        println(v," of type: ", typeof(v)," becomes: ", val, " of type: ", typeof(val))
        
        v = "[1.2, 3.0, 5.98]"
        val = parseTuple(v)
        println(v," of type: ", typeof(v)," becomes: ", val, " of type: ", typeof(val))

        v = "[('mc', 'bc', [1.2, 3.5, 40.9])]"
        val = parseTuple(v)
        println(v," of type: ", typeof(v)," becomes: ", val, " of type: ", typeof(val))
        
    end
    
    function test()
        fn = tempname()
        cfg = """
key: "value"
key2:              ##comment
                   ##comment
    key21: "value" ## comment
                   ##comment
    key22: [1,2,3]
    key23: 234  #comment
    """
        open(fn, "w") do tmp_filehandle
            write(tmp_filehandle, cfg)
        end
        cfg = split(cfg, "\n")

        println("=== Test:  Original configuration ===")
        printlines(cfg)
        # println(cf)
        println("============\n")
        data = readConfigFile(fn)
        println("=== Test: Dict Result ===")
        println(data)
        println("============\n")
        rm(fn)
    end
    
    function pgtest()
        fn = tempname()
        cfg = """
key: "value"
key2:              ##comment
                   ##comment
    key21: "value" ## comment
                   ##comment
    key22: [1,2,3]
    key23: 234  #comment
    """
        open(fn, "w") do tmp_filehandle
            write(tmp_filehandle, cfg)
        end
        cfg = split(cfg, "\n")
        
        println("=== Test:  Original configuration ===")
        printlines(cfg)
        # println(cf)
        println("============\n")

        println("============\n")
        data = pgc.readConfigFile(fn)
        println("=== Test: PG Result ===")
        println(data)
        println("============\n")
        
        # # now do the julia version
        # ln, dataj = readConfigFile(fn)
        # println("=== Test: Julia Result ===")
        # println(dataj)
        # println("============\n")
        #
        # println(data == dataj)
        
        rm(fn)
    end
    
    function printdict(d)
        for (k, v) in d
            if typeof(v) == Dict{Any, Any}
                printdict(v)
            else
                println("key: ", k, "  values: ", v)
            end
        end
    end
    function test_index()
        fn = "/Users/pbmanis/Desktop/2018.09.27_000/ImageSequence_000/.index"
        fn = "/Users/pbmanis/Desktop/Python/mrk-nf107/data_for_testing/CCIV/.index"
        ln, data = readConfigFile(fn)
        printdict(data)
    end
 
end



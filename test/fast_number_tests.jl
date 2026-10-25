# Tests for the tokenizer and number parsing in parser.jl. On well-formed input
# every path gives the same parsed output, so these test the pieces directly:
#
#   1. `_fast_float` and `_fast_int` are bit-identical to `Base.parse`.
#   2. `iter_num` returns the same value and next position as `iter_ws` followed
#      by `fast_parse`. A wrong position would shift every later field of a row.
#   3. `next_line` splits a file into the same lines as the line-vector parser did.
using ExaPowerIO, Test, Random
using ExaPowerIO: _fast_float, _fast_int, fast_parse, iter_num, iter_ws, WordedString

sub(s::String) = SubString(s, 1, ncodeunits(s))

@testset "fast float is bit-identical to Base.parse" begin
    # `===` compares bit patterns, so -0.0 and 0.0 differ.
    edges = ["0", "-0", "0.0", "-0.0", "1", "-1", "1.", ".5", "-.5", "+1.5",
             "1e0", "1E0", "1e-3", "1e+3", "-2.5e-4", "1e22", "1e23", "1e-22",
             "1e-23", "9007199254740992", "9007199254740993", "5e-324",
             "1.7976931348623157e308", "1e400", "1e-400", "123456789012345678",
             "1234567890123456789", "0.1", "0.2", "0.3",
             "00000000000000000000.5", "0.00000000000000000001"]
    for s in edges
        v = _fast_float(sub(s))
        v === nothing && continue          # declined -> Base.parse is used
        @test v === parse(Float64, s)
    end
    # forms the fast path must leave to Base.parse
    for s in ["", "-", "+", ".", "e5", "1e", "1e+", "0x1p3", "Inf", "-Inf",
              "NaN", "nan", "1.5f0", "1 ", " 1", "1;", "--1", "1/3"]
        @test _fast_float(sub(s)) === nothing
    end
    rng = MersenneTwister(20260821)
    for _ in 1:20_000
        ip = rand(rng, 0:10^rand(rng, 1:9))
        fn = rand(rng, 0:9)
        s = fn == 0 ? string(ip) : string(ip, ".", lpad(rand(rng, 0:10^fn-1), fn, '0'))
        rand(rng, Bool) && (s = "-" * s)
        rand(rng) < 0.25 && (s *= "e" * (rand(rng, Bool) ? "-" : "") * string(rand(rng, 0:30)))
        v = _fast_float(sub(s))
        v === nothing || @test v === parse(Float64, s)
    end
    # random bit patterns, which mostly fall outside the fast path
    for _ in 1:5_000
        x = reinterpret(Float64, rand(rng, UInt64))
        isfinite(x) || continue
        v = _fast_float(sub(string(x)))
        v === nothing || @test v === parse(Float64, string(x))
    end
end

@testset "fast int is identical to Base.parse" begin
    rng = MersenneTwister(0xC0FFEE)
    for _ in 1:20_000
        s = string(rand(rng, Int64) ÷ rand(rng, 1:10^6))
        v = _fast_int(sub(s))
        v === nothing || @test v == parse(Int, s)
    end
    for s in ["", "-", "+", "1.0", "1e3", "0x10", "1 ", " 1", "--1",
              "1234567890123456789012"]
        @test _fast_int(sub(s)) === nothing
    end
end

@testset "fused iter_num matches iter_ws + fast_parse (value and position)" begin
    function check(line::AbstractString)
        l = sub(String(line))
        ws = WordedString(l, l.ncodeunits)
        pos = 1
        while true
            tok, nxt = iter_ws(ws, pos)
            ref = tryparse(Float64, tok)
            ref === nothing && break        # ';', ']': both paths throw
            got = iter_num(Float64, ws, pos)
            @test got[1] === ref
            @test got[2] == nxt
            nxt == 0 && break
            pos = nxt
        end
    end
    for f in readdir(joinpath(@__DIR__, "..", "data"))
        f == "LICENSE" && continue
        for l in eachline(joinpath(@__DIR__, "..", "data", f))
            isempty(l) || check(l)
        end
    end
    for l in ["1 2 3", "  1.5   -2.5e-3  ;", "1;", "-0.0 0.0", "1e 2", "--1 3",
              "12345678901234567890.5 1", "1e400 2", "1e-400 2", "0x10 1",
              "Inf NaN 1", "1.5\t2.5\t3", "1", "1 ", " 1", "1]2",
              "9007199254740993 1", "1e22 1e23", "1.5\u00a02", "1é 2"]
        check(l)
    end
end

@testset "a row with too few fields raises" begin
    # The row macro calls iter_num with start == 0 when a row has fewer fields
    # than its section declares; this raises the error parse gives on "".
    ws = WordedString(sub("1"), 1)
    @test_throws ArgumentError iter_num(Float64, ws, 0)
    @test_throws ArgumentError iter_num(Int, ws, 0)
end

@testset "only ASCII bytes are delimiters" begin
    # 0x85 and 0xa0 are `isspace` as Chars but continuation bytes in UTF-8
    @test !any(ExaPowerIO.SPACE_B[129:end])
    @test !any(ExaPowerIO.ENDS_B[129:end])
end

@testset "next_line matches split(replace(s, \"\\r\\n\" => \"\\n\"), \"\\n\")" begin
    for s in ["", "\n", "\n\n", "a", "a\n", "a\nb", "a\r\nb\r\n", "\r\n", "a\r",
              "a\r\r\n", "x é\ny é", "é\r\nü", "% ü"]
        lines = SubString{String}[]
        pos = 1
        while (nxt = ExaPowerIO.next_line(s, pos, ncodeunits(s))) !== nothing
            line, pos = nxt
            push!(lines, line)
        end
        @test lines == split(replace(s, "\r\n" => "\n"), "\n")
    end
end

@testset "CRLF and non-ASCII comments parse as the original file" begin
    src = joinpath(@__DIR__, "..", "data", "pglib_opf_case5_pjm_storage.m")
    text = replace(read(src, String), "\r\n" => "\n")
    variant = replace(join((startswith(l, "%") ? l * " é" : l for l in split(text, "\n")), "\n"),
                      "\n" => "\r\n")
    path, io = mktemp()
    write(io, variant); close(io)
    a = ExaPowerIO.parse_matpower(src)
    b = ExaPowerIO.parse_matpower(path)
    for f in fieldnames(ExaPowerIO.PowerData)
        @test getfield(a, f) == getfield(b, f)
    end
    rm(path)
end

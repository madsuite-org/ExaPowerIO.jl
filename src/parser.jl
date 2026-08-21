# FOR CONVENIENCE, the MATPOWER file spec
# https://matpower.app/manual/matpower/DataFileFormat.html

"""
    struct BusData{T <: Real}
        i :: Int
        bus_i :: Int
        type :: Int
        pd :: T
        qd :: T
        gs :: T
        bs :: T
        area :: Int
        vm :: T
        va :: T
        baseKV :: T
        zone :: Int
        vmax :: T
        vmin :: T
    end
"""
struct BusData{T <: Real}
    i :: Int
    bus_i :: Int
    type :: Int
    pd :: T
    qd :: T
    gs :: T
    bs :: T
    area :: Int
    vm :: T
    va :: T
    baseKV :: T
    zone :: Int
    vmax :: T
    vmin :: T
end

"""
    struct BranchData{T <: Real}
        i :: Int
        f_bus :: Int
        t_bus :: Int
        br_r :: T
        br_x :: T
        b_fr :: T,
        b_to :: T,
        g_fr :: T,
        g_to :: T,
        rate_a ::T
        rate_b :: T
        rate_c :: T
        tap :: T
        shift :: T
        status :: Int
        angmin :: T
        angmax :: T
        f_idx::Int
        t_idx::Int
        c1 :: T
        c2 :: T
        c3 :: T
        c4 :: T
        c5 :: T
        c6 :: T
        c7 :: T
        c8 :: T
    end

f_bus and t_bus are indices into the PowerData.bus Vector, not bus_i values
"""
struct BranchData{T <: Real}
    i :: Int
    f_bus :: Int
    t_bus :: Int
    br_r :: T
    br_x :: T
    b_fr :: T
    b_to :: T
    g_fr :: T
    g_to :: T
    rate_a ::T
    rate_b :: T
    rate_c :: T
    tap :: T
    shift :: T
    status :: Int
    angmin :: T
    angmax :: T
    f_idx::Int
    t_idx::Int
    c1 :: T
    c2 :: T
    c3 :: T
    c4 :: T
    c5 :: T
    c6 :: T
    c7 :: T
    c8 :: T
end
function BranchData{T}(
    i::Int,
    f_bus::Int,
    t_bus::Int,
    br_r::T,
    br_x::T,
    b_fr::T,
    b_to::T,
    g_fr::T,
    g_to::T,
    rate_a::T,
    rate_b::T,
    rate_c::T,
    tap::T,
    shift::T,
    status::Int,
    angmin::T,
    angmax::T,
    f_idx::Int,
    t_idx::Int
) where {T<:Real}
    x = br_r + im * br_x
    xi = inv(x)
    y = ifelse(isfinite(xi), xi, zero(xi))
    g = real(y)
    b = imag(y)
    if isapprox(tap, T(0.0))
        tap = T(1.0)
    end
    tr = tap * cos(shift)
    ti = tap * sin(shift)
    ttm = tr^2 + ti^2
    c1 = (-g * tr - b * ti) / ttm
    c2 = (-b * tr + g * ti) / ttm
    c3 = (-g * tr + b * ti) / ttm
    c4 = (-b * tr - g * ti) / ttm
    c5 = (g + g_fr) / ttm
    c6 = (b + b_fr) / ttm
    c7 = (g + g_to)
    c8 = (b + b_to)
    BranchData{T}(
        i,
        f_bus,
        t_bus,
        br_r,
        br_x,
        b_fr,
        b_to,
        g_fr,
        g_to,
        rate_a,
        rate_b,
        rate_c,
        tap,
        shift,
        status,
        angmin,
        angmax,
        f_idx,
        t_idx,
        c1,
        c2,
        c3,
        c4,
        c5,
        c6,
        c7,
        c8,
    )
end

"""
    struct ArcData{T <: Real}
        i :: Int
        bus :: Int
        rate_a :: T
    end
"""
struct ArcData{T <: Real}
    i :: Int
    bus :: Int
    rate_a :: T
end

"""
    struct StorageData{T <: Real}
        i :: Int
        storage_bus :: Int
        ps :: T
        qs :: T
        energy :: T
        energy_rating :: T
        charge_rating :: T
        discharge_rating :: T
        charge_efficiency :: T
        discharge_efficiency :: T
        thermal_rating :: T
        qmin :: T
        qmax :: T
        r :: T
        x :: T
        p_loss :: T
        q_loss :: T
        status :: Int
    end
"""
struct StorageData{T <: Real}
    i :: Int
    storage_bus :: Int
    Pexts :: T
    Qexts :: T
    energy :: T
    energy_rating :: T
    charge_rating :: T
    discharge_rating :: T
    charge_efficiency :: T
    discharge_efficiency :: T
    thermal_rating :: T
    qmin :: T
    qmax :: T
    Zr :: T
    Zim :: T
    p_loss :: T
    q_loss :: T
    status :: Int
end

"""
    struct GenData{T <: Real}
        i :: Int
        bus :: Int
        pg :: T
        qg :: T
        qmax :: T
        qmin :: T
        vg :: T
        mbase :: T
        status :: Int
        pmax :: T
        pmin :: T
        model_poly :: Bool
        startup :: T
        shutdown :: T
        n :: Int
        c :: NTuple{3, T}
    end

bus is an index into the PowerData.bus Vector, not bus_i values
"""
struct GenData{T <: Real}
    i :: Int
    bus :: Int
    pg :: T
    qg :: T
    qmax :: T
    qmin :: T
    vg :: T
    mbase :: T
    status :: Int
    pmax :: T
    pmin :: T
    model_poly :: Bool
    startup :: T
    shutdown :: T
    n :: Int
    c :: NTuple{3, T}
end

"""
    struct Data{T <: Real}
        version :: String
        baseMVA :: T
        bus :: Vector{BusData{T}}
        gen :: Vector{GenData{T}}
        branch :: Vector{BranchData{T}}
        storage :: Vector{StorageData{T}}
    end
```version```, ```baseMVA```, ```bus```, ```gen```, ```branch```, and ```storage```
all corespond to members of the mpc object created by a matpower file.
Their fields correspond exactly with the columns of the relevant ```mpc``` member.
"""
struct PowerData{
    T <: Real,
    VBusT <: AbstractVector{BusData{T}},
    VGenT <: AbstractVector{GenData{T}},
    VBranchT <: AbstractVector{BranchData{T}},
    VStorageT <: AbstractVector{StorageData{T}},
    VArcT <: AbstractVector{ArcData{T}}
}
    version :: String
    baseMVA :: T
    bus :: VBusT
    gen :: VGenT
    branch :: VBranchT
    arc :: VArcT
    storage :: VStorageT
end

const MATPOWER_KEYS :: Vector{String} = ["version", "baseMVA", "areas", "bus", "gencost", "gen", "branch", "storage"];
const NULL_VIEW::SubString{String} = SubString("", 1, 0)
is_end(c::Char) = isspace(c) || c in "=;[]%"
const KEY_MIN_LEN = 4

struct WordedString
    s :: SubString{String}
    len :: Int
end

# Byte-indexed delimiter tables, indexed by byte+1 so byte 0x00 stays in range.
# These beat the equivalent chain of byte comparisons by 7.9 ns/token measured
# on case9241 (`bench/tokprobe.jl`, variant A vs B) -- the table is one load,
# the comparison chain is a branch per delimiter. Keep the table.
const ENDS_B = ntuple(i -> is_end(Char(i - 1)), 256)
const SPACE_B = ntuple(i -> isspace(Char(i - 1)), 256)

# Tokenizer over raw code units. The previous version indexed the SubString by
# CHARACTER (`ws.s[left]`), decoding UTF-8 on every access, and narrowed each
# char through `Int8` -- which also throws on any byte > 127. MATPOWER numeric
# rows are ASCII, so byte indexing is equivalent here and materially cheaper.
# The bounds test now precedes the read: the old `isspace(ws.s[left]) && left
# <= ws.len` read first and checked after, which under @inbounds reads one past
# the end on a line of trailing blanks.
@inline @inbounds function iter_ws(ws :: WordedString, start :: Int) :: Tuple{SubString{String}, Int}
    len = ws.len
    if start > len
        return (NULL_VIEW, 0)
    end
    str = ws.s
    left = start
    while left <= len && SPACE_B[codeunit(str, left) + 1]
        left += 1
    end
    if left > len || codeunit(str, left) == 0x25   # '%'
        return (NULL_VIEW, 0)
    end
    right = left
    while right <= len && !ENDS_B[codeunit(str, right) + 1]
        right += 1
    end
    # right is non-inclusive
    if ENDS_B[codeunit(str, left) + 1]
        right += 1
    end
    (SubString(str, left, right - 1), right)
end

# --- exact fast number parsing ----------------------------------------------
# Clinger's fast path: when the significand fits in 2^53 and the decimal
# exponent is within +/-22, BOTH the significand and 10^|e| are exactly
# representable in Float64, so one IEEE multiply (or divide) is correctly
# rounded -- bit-identical to Base.parse, not an approximation. Every other
# shape -- oversized mantissa, big exponent, hex, Inf/NaN, trailing garbage --
# returns `nothing` and falls back to Base.parse. So this path can only be
# faster, never different; `bench/fuzz_numbers.jl` checks that claim directly.
const POW10_EXACT = ntuple(i -> 10.0^(i - 1), 23)   # 10^0 .. 10^22, each exact

@inline function _fast_float(s::SubString{String})
    n = ncodeunits(s)
    n == 0 && return nothing
    i = 1
    b = codeunit(s, i)
    neg = false
    if b == 0x2d          # '-'
        neg = true; i += 1
    elseif b == 0x2b      # '+'
        i += 1
    end
    i > n && return nothing
    mant = UInt64(0); ndig = 0; frac = 0; seen = false
    while i <= n
        b = codeunit(s, i)
        (b < 0x30 || b > 0x39) && break
        ndig += 1
        ndig > 18 && return nothing      # bail out before UInt64 could overflow
        mant = mant * 10 + (b - 0x30)
        seen = true; i += 1
    end
    if i <= n && codeunit(s, i) == 0x2e  # '.'
        i += 1
        while i <= n
            b = codeunit(s, i)
            (b < 0x30 || b > 0x39) && break
            ndig += 1
            ndig > 18 && return nothing
            mant = mant * 10 + (b - 0x30)
            frac += 1; seen = true; i += 1
        end
    end
    seen || return nothing
    exp10 = -frac
    if i <= n
        b = codeunit(s, i)
        if b == 0x65 || b == 0x45        # 'e' / 'E'
            i += 1
            i > n && return nothing
            b = codeunit(s, i)
            eneg = false
            if b == 0x2d
                eneg = true; i += 1
            elseif b == 0x2b
                i += 1
            end
            i > n && return nothing
            ev = 0; eseen = false
            while i <= n
                b = codeunit(s, i)
                (b < 0x30 || b > 0x39) && break
                ev = ev * 10 + Int(b - 0x30)
                ev > 1000 && return nothing
                eseen = true; i += 1
            end
            eseen || return nothing
            exp10 += eneg ? -ev : ev
        end
    end
    i <= n && return nothing             # trailing garbage -> let Base decide
    mant > (UInt64(1) << 53) && return nothing
    m = Float64(mant)
    if exp10 == 0
        v = m
    elseif 0 < exp10 <= 22
        v = m * POW10_EXACT[exp10 + 1]
    elseif -22 <= exp10 < 0
        v = m / POW10_EXACT[-exp10 + 1]
    else
        return nothing
    end
    return neg ? -v : v
end

@inline function _fast_int(s::SubString{String})
    n = ncodeunits(s)
    n == 0 && return nothing
    i = 1
    b = codeunit(s, i)
    neg = false
    if b == 0x2d
        neg = true; i += 1
    elseif b == 0x2b
        i += 1
    end
    i > n && return nothing
    v = 0; nd = 0
    while i <= n
        b = codeunit(s, i)
        (b < 0x30 || b > 0x39) && return nothing
        nd += 1
        nd > 18 && return nothing
        v = v * 10 + Int(b - 0x30)
        i += 1
    end
    nd == 0 && return nothing
    return neg ? -v : v
end

# Only Float64 and Int get a fast path; every other T goes straight to
# Base.parse, so generic-precision parsing (Float32, BigFloat, ...) is
# bit-for-bit untouched by this.
# Fused scan-and-parse. `iter_ws` walks a token's bytes to find where it ends,
# and `_fast_float` then walks the SAME bytes again to read the digits. For the
# numeric rows -- which are ~93% of all tokens in a MATPOWER file -- one pass
# can do both: accumulate the significand while scanning, and stop at the
# delimiter that would have ended the token anyway.
#
# `nextpos` is exactly what `iter_ws` would have returned, so callers cannot
# tell the two apart. Any shape the fused path does not handle (comment, empty,
# oversized mantissa, exponent out of range, trailing garbage) rewinds and
# takes the original `iter_ws` + `fast_parse` route, so behaviour is identical
# by construction -- the fused path is an accelerator, never a second parser
# with its own opinion. `bench/fuse_equiv.jl` checks value AND nextpos agree
# token-for-token over every pglib case.
# `start == 0` means the previous field consumed the last token on the line.
# The old macro special-cased that with `state[2] == 0 ? (NULL_VIEW, 0) : ...`,
# so a short row raised "cannot parse \"\" as Float64" from Base. Reproduce that
# exactly: without this guard we would index the line at offset 0 under
# @inbounds, which reads out of bounds instead of raising. A well-formed file
# never reaches here, so no amount of pglib data would have caught it.
@inline function iter_num(::Type{T}, ws::WordedString, start::Int) where {T}
    start < 1 && return (fast_parse(T, NULL_VIEW), 0)
    tok, nxt = iter_ws(ws, start)
    return (fast_parse(T, tok), nxt)
end

@inline function iter_num(::Type{Float64}, ws::WordedString, start::Int)
    @inbounds begin
        len = ws.len
        (start < 1 || start > len) && return _iter_num_slow(Float64, ws, start)
        str = ws.s
        i = start
        while i <= len && SPACE_B[codeunit(str, i) + 1]
            i += 1
        end
        (i > len || codeunit(str, i) == 0x25) && return _iter_num_slow(Float64, ws, start)
        b = codeunit(str, i)
        neg = false
        if b == 0x2d
            neg = true; i += 1
        elseif b == 0x2b
            i += 1
        end
        i > len && return _iter_num_slow(Float64, ws, start)
        mant = UInt64(0); ndig = 0; frac = 0; seen = false
        while i <= len
            b = codeunit(str, i)
            (b < 0x30 || b > 0x39) && break
            ndig += 1
            ndig > 18 && return _iter_num_slow(Float64, ws, start)
            mant = mant * 10 + (b - 0x30)
            seen = true; i += 1
        end
        if i <= len && codeunit(str, i) == 0x2e
            i += 1
            while i <= len
                b = codeunit(str, i)
                (b < 0x30 || b > 0x39) && break
                ndig += 1
                ndig > 18 && return _iter_num_slow(Float64, ws, start)
                mant = mant * 10 + (b - 0x30)
                frac += 1; seen = true; i += 1
            end
        end
        seen || return _iter_num_slow(Float64, ws, start)
        exp10 = -frac
        if i <= len
            b = codeunit(str, i)
            if b == 0x65 || b == 0x45
                i += 1
                i > len && return _iter_num_slow(Float64, ws, start)
                b = codeunit(str, i)
                eneg = false
                if b == 0x2d
                    eneg = true; i += 1
                elseif b == 0x2b
                    i += 1
                end
                i > len && return _iter_num_slow(Float64, ws, start)
                ev = 0; eseen = false
                while i <= len
                    b = codeunit(str, i)
                    (b < 0x30 || b > 0x39) && break
                    ev = ev * 10 + Int(b - 0x30)
                    ev > 1000 && return _iter_num_slow(Float64, ws, start)
                    eseen = true; i += 1
                end
                eseen || return _iter_num_slow(Float64, ws, start)
                exp10 += eneg ? -ev : ev
            end
        end
        # the token must end HERE, at a real delimiter (or end of line), or the
        # original tokenizer would have taken more bytes than we just consumed
        if i <= len && !ENDS_B[codeunit(str, i) + 1]
            return _iter_num_slow(Float64, ws, start)
        end
        mant > (UInt64(1) << 53) && return _iter_num_slow(Float64, ws, start)
        m = Float64(mant)
        if exp10 == 0
            v = m
        elseif 0 < exp10 <= 22
            v = m * POW10_EXACT[exp10 + 1]
        elseif -22 <= exp10 < 0
            v = m / POW10_EXACT[-exp10 + 1]
        else
            return _iter_num_slow(Float64, ws, start)
        end
        return (neg ? -v : v, i)
    end
end

@noinline function _iter_num_slow(::Type{T}, ws::WordedString, start::Int) where {T}
    start < 1 && return (fast_parse(T, NULL_VIEW), 0)
    tok, nxt = iter_ws(ws, start)
    return (fast_parse(T, tok), nxt)
end

@inline fast_parse(::Type{T}, s::SubString{String}) where {T} = parse(T, s)
@inline function fast_parse(::Type{Float64}, s::SubString{String})
    v = _fast_float(s)
    return v === nothing ? parse(Float64, s) : v
end
@inline function fast_parse(::Type{Int}, s::SubString{String})
    v = _fast_int(s)
    return v === nothing ? parse(Int, s) : v
end

macro iter_to_ntuple(n, iter_expr, types)
    iter_sym = gensym("iter")
    state_sym = gensym("state")
    x_syms = [gensym("x") for _ in 1:n]

    body = Expr[]
    push!(body, :($iter_sym = $(esc(iter_expr))))
    pos_sym = gensym("pos")
    push!(body, :($pos_sym = 1))

    length(types.args) != n && error("types provided to @iter_to_ntuple had length $(length(types.args)) instead of $n")
    for i in 1:n
        # one call per field: scans the token and parses it in a single pass
        push!(body, :($state_sym = iter_num($(esc(types.args[i])), $iter_sym, $pos_sym)))
        push!(body, :($(x_syms[i]) = $state_sym[1]))
        if i < n
            push!(body, :($pos_sym = $state_sym[2]))
        end
    end
    push!(body, Expr(:tuple, x_syms...))

    return Expr(:block, body...)
end

# Walk the '\n'-separated lines of `s` without materializing them. Semantics
# match `split(s, "\n")` exactly -- including the single empty trailing field
# when `s` ends in a newline -- so line numbering and array lengths are
# unchanged. `\r` is left on the line, as `split` leaves it, so CRLF files
# behave exactly as before.
@inline function next_line(s::String, pos::Int, n::Int)
    if pos > n
        # `split` emits one empty final field iff the input ended with a
        # newline (and one for the empty input); anything past that is done.
        (pos == n + 1 && (n == 0 || codeunit(s, n) == 0x0a)) &&
            return (SubString(s, n + 1, n), n + 2)
        return nothing
    end
    # Search the CODE UNITS, not the String. `findnext('\n', ::String, ...)`
    # goes through the character-oriented path; the byte form dispatches to
    # memchr and is 29% faster per line (41.1 vs 57.7 ns, bench/lineprobe.jl).
    # A hand-rolled byte loop is much SLOWER (92 ns) -- it gives up the SIMD
    # memchr -- so this is the fast form, not an intermediate step.
    nl = findnext(==(0x0a), codeunits(s), pos)
    nl === nothing && return (SubString(s, pos, n), n + 1)
    return (SubString(s, pos, nl - 1), nl + 1)
end

# `pos` is the byte offset of the line FOLLOWING the "mpc.<key> = [" line,
# which is where the original counted from. Returns the number of data rows
# before the closing "]". Indexing an empty line raises, exactly as before.
function get_arr_len(s::String, n::Int, pos::Int)::Int
    count = 0
    p = pos
    while true
        nxt = next_line(s, p, n)
        nxt === nothing && error("Array beginning at byte $pos was not closed")
        line, p = nxt
        line[1] == ']' && return count
        count += 1
    end
end

@inbounds @inline @views function parse_matpower_inner(::Type{T}, ::Type{V}, fname :: String, filtered :: Bool) where {T<:Real, V<:AbstractVector}
    # `read(open(fname), String)` never closed the stream; `read(fname, String)`
    # does. The line vector is gone -- the parser only ever walks forward.
    # Measured on case9241 (4.75 MiB, node 000): dropping it saves 5.03 MiB of
    # the 15.54 MiB the parse allocated, and 3.1 ms of 89.6 ms. The memory is
    # the point here; the time is a rounding error next to tokenizing.
    fstring = read(fname, String)
    nbytes = ncodeunits(fstring)
    in_array = false
    cur_key = ""
    bus = BusData{T}[]
    gen = GenData{T}[]
    skipped_gens = Int[]
    branch = BranchData{T}[]
    storage = StorageData{T}[]

    row_num = 1
    version = ""
    baseMVA :: T = T(0.0)
    # `Int[]`, not `[]`. An untyped literal is a `Vector{Any}`, and the
    # annotation then converts it — `convert(Vector{Int64}, ::Vector{Any})`
    # copies element by element through a dynamic `setindex!`, which
    # `juliac --trim=safe` refuses. Building it at the right element type
    # costs nothing and removes the conversion entirely.
    bus_map :: Vector{Int} = Int[]
    bus_offset :: Int = 0
    line_ind = 0
    num_branch = 0
    cur_branch = 1
    biggest_gen = -1
    biggest_gen_pmax = -Inf
    num_skipped_gens = 0
    cur_skipped_gen = 1
    pos = 1
    while true
        nxt = next_line(fstring, pos, nbytes)
        nxt === nothing && break
        # advance the cursor BEFORE the body -- the body uses `continue`, and
        # `pos` must already point at the next line for `get_arr_len` anyway.
        line, pos = nxt
        line_len = line.ncodeunits
        line_ind += 1
        line_len != 0 && codeunit(line, 1) == 0x25 && continue   # '%'
        if in_array && line_len != 0
            # exactly `startswith(line, "];")`, as two byte compares -- this runs
            # once per data row, not once per array
            if line_len >= 2 && codeunit(line, 1) == 0x5d && codeunit(line, 2) == 0x3b
                if cur_key == "bus"
                    bus_offset = minimum(b -> b.bus_i, bus) - 1
                    max_bus = maximum(b -> b.bus_i, bus)
                    bus_map = [0 for _ in 1:(max_bus-bus_offset)]
                    for (i, b) in enumerate(bus)
                        bus_map[b.bus_i - bus_offset] = i
                    end
                end
                in_array = false
            elseif cur_key == "bus"
                bus_words = @iter_to_ntuple 13 WordedString(line, line_len) (Int, Int, T, T, T, T, Int, T, T, T, Int, T, T)
                bus[row_num] = BusData(
                    row_num,
                    bus_words[1],
                    bus_words[2],
                    bus_words[3] / baseMVA,
                    bus_words[4] / baseMVA,
                    bus_words[5] / baseMVA,
                    bus_words[6] / baseMVA,
                    bus_words[7],
                    bus_words[8],
                    bus_words[9],
                    bus_words[10],
                    bus_words[11],
                    bus_words[12],
                    bus_words[13],
                )
                if filtered && bus[row_num].type == 4
                    pop!(bus)
                    continue
                end
            elseif cur_key == "gen"
                gen_words = @iter_to_ntuple 10 WordedString(line, line_len) (Int, T, T, T, T, T, T, Int, T, T)
                if gen_words[8] != 0 && gen_words[10] > biggest_gen_pmax
                    biggest_gen_pmax = gen_words[10]
                    biggest_gen = row_num
                end
                gen[row_num] = GenData(
                    row_num,
                    bus_map[gen_words[1] - bus_offset],
                    gen_words[2] / baseMVA,
                    gen_words[3] / baseMVA,
                    gen_words[4] / baseMVA,
                    gen_words[5] / baseMVA,
                    gen_words[6],
                    gen_words[7],
                    gen_words[8],
                    gen_words[9] / baseMVA,
                    gen_words[10] / baseMVA,
                    false,
                    T(0),
                    T(0),
                    0,
                    (T(0), T(0), T(0)),
                )
                if filtered && gen[row_num].status == 0
                    pop!(gen)
                    num_skipped_gens += 1
                    skipped_gens[num_skipped_gens] = row_num
                    skipped_gens[num_skipped_gens+1] = 0
                    continue
                end
            elseif cur_key == "gencost"
                genc_words = @iter_to_ntuple 7 WordedString(line, line_len) (Int, T, T, Int, T, T, T)
                if filtered && row_num == skipped_gens[cur_skipped_gen]
                    cur_skipped_gen += 1
                    continue
                end
                model_poly = genc_words[1] == 2
                n = genc_words[4]
                normalize_cost = let baseMVA = baseMVA
                    function normalize_cost(i :: Int)
                        c = genc_words[4 + i]
                        return model_poly ? baseMVA ^ (n-i) * c : c
                    end
                end
                gen[row_num] = GenData(
                    row_num,
                    gen[row_num].bus,
                    gen[row_num].pg,
                    gen[row_num].qg,
                    gen[row_num].qmax,
                    gen[row_num].qmin,
                    gen[row_num].vg,
                    gen[row_num].mbase,
                    gen[row_num].status,
                    gen[row_num].pmax,
                    gen[row_num].pmin,
                    model_poly,
                    genc_words[2],
                    genc_words[3],
                    n,
                    ntuple(normalize_cost, 3)
                )
            elseif cur_key == "branch"
                branch_words = @iter_to_ntuple 13 WordedString(line, line_len) (Int, Int, T, T, T, T, T, T, T, T, Int, T, T)
                branch[row_num] = BranchData{T}(
                    row_num,
                    bus_map[branch_words[1] - bus_offset],
                    bus_map[branch_words[2] - bus_offset],
                    branch_words[3],
                    branch_words[4],
                    branch_words[5] / T(2.0),
                    branch_words[5] / T(2.0),
                    T(0.0),
                    T(0.0),
                    branch_words[6] / baseMVA,
                    branch_words[7] / baseMVA,
                    branch_words[8] / baseMVA,
                    branch_words[9],
                    (branch_words[10]) / T(180.0) * T(pi),
                    branch_words[11],
                    branch_words[12] / T(180.0) * T(pi),
                    branch_words[13] / T(180.0) * T(pi),
                    cur_branch,
                    cur_branch + num_branch
                )
                if filtered && branch[row_num].status == 0
                    pop!(branch)
                    continue
                end
                cur_branch += 1
            elseif cur_key == "storage"
                storage_words = @iter_to_ntuple 17 WordedString(line, line_len) (Int, T, T, T, T, T, T, T, T, T, T, T, T, T, T, T, Int)
                storage[row_num] = StorageData(
                    row_num,
                    storage_words[1],
                    storage_words[2],
                    storage_words[3],
                    storage_words[4] / baseMVA,
                    storage_words[5] / baseMVA,
                    storage_words[6] / baseMVA,
                    storage_words[7] / baseMVA,
                    storage_words[8],
                    storage_words[9],
                    storage_words[10] / baseMVA,
                    storage_words[11] / baseMVA,
                    storage_words[12] / baseMVA,
                    storage_words[13],
                    storage_words[14],
                    storage_words[15],
                    storage_words[16],
                    storage_words[17],
                )
            end
            if in_array
                row_num += 1
            end
        elseif line_len > 0 && !startswith(line, "function") && isletter(line[1])
            cur_key = ""
            for key in MATPOWER_KEYS
                full_name = "mpc.$key ="
                if startswith(line, full_name)
                    cur_key = key
                    break
                end
            end

            if cur_key == ""
                @warn "Unrecognized assignment at line $line_ind.\n\t$line" 
                continue
            end
            if cur_key == "version"
                raw_data = split(line)[3]
                version = String(raw_data[2:raw_data.ncodeunits-2])
            elseif cur_key == "baseMVA"
                word = split(line)[3]
                baseMVA = parse(T, word[1:length(word)-1]) :: T
            else
                arr_len = get_arr_len(fstring, nbytes, pos)
                if cur_key == "bus"
                    bus = V{BusData{T}}(undef, arr_len)
                elseif cur_key == "gen"
                    gen = V{GenData{T}}(undef, arr_len)
                    if filtered
                        # extra elem for terminator
                        skipped_gens = Vector{Int}(undef, arr_len + 1)
                        skipped_gens[1] = 0
                    end
                elseif cur_key == "branch"
                    branch = V{BranchData{T}}(undef, arr_len)
                    num_branch = arr_len
                elseif cur_key == "storage"
                    storage = V{StorageData{T}}(undef, arr_len)
                end
                in_array = true
                row_num = 1
            end
        end
    end

    has_gen = [false for _ in 1:length(bus)]
    look_for_ref = false
    for gen in gen
        if gen.status == 1
            has_gen[gen.bus] = true
        end
    end
    for (i, b) in enumerate(bus)
        if has_gen[i] && b.type == 1
            bus[i] = BusData(
                bus[i].i,
                b.bus_i,
                2,
                b.pd,
                b.qd,
                b.gs,
                b.bs,
                b.area,
                b.vm,
                b.va,
                b.baseKV,
                b.zone,
                b.vmax,
                b.vmin
            )
        elseif !has_gen[i] && (b.type == 2 || b.type == 3)
            if bus[i].type == 3
                look_for_ref = true
            end
            bus[i] = BusData(
                bus[i].i,
                b.bus_i,
                1,
                b.pd,
                b.qd,
                b.gs,
                b.bs,
                b.area,
                b.vm,
                b.va,
                b.baseKV,
                b.zone,
                b.vmax,
                b.vmin
            )
        end
    end
    if look_for_ref
        b = bus[gen[biggest_gen].bus]
        bus[gen[biggest_gen].bus] = BusData(
            b.i,
            b.bus_i,
            3,
            b.pd,
            b.qd,
            b.gs,
            b.bs,
            b.area,
            b.vm,
            b.va,
            b.baseKV,
            b.zone,
            b.vmax,
            b.vmin
        )
    end

    # `num_branch` still holds the DECLARED row count here; `t_idx` was computed
    # against it as `cur_branch + num_branch` while rows were read. Filtering
    # `pop!`s inactive branches, so if any were dropped the final count differs
    # and every `t_idx` is stale -- that, and only that, is why the rewrite
    # below exists. When nothing was dropped the values are already right, and
    # the rewrite is pure cost: it rebuilds a 216-byte struct per branch AND
    # recomputes all eight admittance coefficients through the 19-argument
    # constructor. Skip it in that case.
    branch_declared = num_branch
    num_branch = length(branch)
    arc = V{ArcData{T}}(undef, num_branch * 2)
    needs_tidx_fix = filtered && num_branch != branch_declared
    # `b` from the iterator already IS branch[i]; the previous version wrote
    # `branch[i].<field>` nineteen times per row, reloading a 27-field struct
    # from the array on every single field access. `filtered` is loop-invariant
    # so it is hoisted out rather than retested once per branch.
    if needs_tidx_fix
        for (i, b) in enumerate(branch)
            branch[i] = BranchData{T}(
                b.i, b.f_bus, b.t_bus, b.br_r, b.br_x, b.b_fr, b.b_to,
                b.g_fr, b.g_to, b.rate_a, b.rate_b, b.rate_c, b.tap, b.shift,
                b.status, b.angmin, b.angmax, b.f_idx, b.f_idx + num_branch,
            )
            arc[i] = ArcData(i, b.f_bus, b.rate_a)
            arc[i+num_branch] = ArcData(i+num_branch, b.t_bus, b.rate_a)
        end
    else
        for (i, b) in enumerate(branch)
            arc[i] = ArcData(i, b.f_bus, b.rate_a)
            arc[i+num_branch] = ArcData(i+num_branch, b.t_bus, b.rate_a)
        end
    end

    return PowerData(version, baseMVA, bus, gen, branch, arc, storage)
end

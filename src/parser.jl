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

# Byte classes for the tokenizer, indexed by `byte + 1`. Only ASCII bytes are
# classified, so a byte inside a multibyte character is never a delimiter
# (`isspace(Char(0x85))` and `isspace(Char(0xa0))` are true, but 0x85 and 0xa0
# are continuation bytes in UTF-8).
const ENDS_B = ntuple(i -> i <= 128 && is_end(Char(i - 1)), 256)
const SPACE_B = ntuple(i -> i <= 128 && isspace(Char(i - 1)), 256)

# The token starting at or after byte `start` of `ws`, and the byte after it.
# Returns `(NULL_VIEW, 0)` at the end of the line or at a `%` comment.
@inline function iter_ws(ws :: WordedString, start :: Int) :: Tuple{SubString{String}, Int}
    len = ws.len
    if start > len
        return (NULL_VIEW, 0)
    end
    str = ws.s
    left = start
    @inbounds while left <= len && SPACE_B[codeunit(str, left) + 1]
        left += 1
    end
    if left > len || @inbounds(codeunit(str, left)) == UInt8('%')
        return (NULL_VIEW, 0)
    end
    right = left
    @inbounds while right <= len && !ENDS_B[codeunit(str, right) + 1]
        right += 1
    end
    # right is non-inclusive
    if @inbounds ENDS_B[codeunit(str, left) + 1]
        right += 1
    end
    (SubString(str, left, thisind(str, right - 1)), right)
end

const POW10_EXACT = ntuple(i -> 10.0^(i - 1), 23)   # 10^0 .. 10^22, each exact

@inline isdigit_b(b::UInt8) = 0x30 <= b <= 0x39

# Reads a decimal number from code units `i:n` of `s` and returns `(value, next)`,
# where `next` is the index of the first byte not read, or `nothing` if no number
# starts at `i` or it cannot be read exactly. It reads exactly when the significand
# has at most 18 digits and is at most 2^53, and the decimal exponent is within ±22.
# Then the significand and 10^|e| are both exact Float64 values, so one multiply or
# divide is correctly rounded and the value is bit-identical to `parse(Float64, ...)`
# (Clinger's fast path). Requires `1 <= i` and `n <= ncodeunits(s)`.
@inline function scan_float(s::SubString{String}, i::Int, n::Int)
    @inbounds begin
        i > n && return nothing
        b = codeunit(s, i)
        neg = b == UInt8('-')
        (neg || b == UInt8('+')) && (i += 1)
        mant = UInt64(0); ndig = 0; frac = 0
        while i <= n && isdigit_b(codeunit(s, i))
            ndig += 1
            ndig > 18 && return nothing   # before mant can overflow UInt64
            mant = mant * 10 + (codeunit(s, i) - 0x30)
            i += 1
        end
        if i <= n && codeunit(s, i) == UInt8('.')
            i += 1
            while i <= n && isdigit_b(codeunit(s, i))
                ndig += 1
                ndig > 18 && return nothing
                mant = mant * 10 + (codeunit(s, i) - 0x30)
                frac += 1; i += 1
            end
        end
        ndig == 0 && return nothing
        exp10 = -frac
        if i <= n && (codeunit(s, i) == UInt8('e') || codeunit(s, i) == UInt8('E'))
            i += 1
            i > n && return nothing
            b = codeunit(s, i)
            eneg = b == UInt8('-')
            (eneg || b == UInt8('+')) && (i += 1)
            ev = 0; edig = 0
            while i <= n && isdigit_b(codeunit(s, i))
                ev = ev * 10 + Int(codeunit(s, i) - 0x30)
                ev > 1000 && return nothing
                edig += 1; i += 1
            end
            edig == 0 && return nothing
            exp10 += eneg ? -ev : ev
        end
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
        return (neg ? -v : v, i)
    end
end

# `s` as a Float64 if `scan_float` reads all of it, else `nothing`.
@inline function _fast_float(s::SubString{String})
    n = ncodeunits(s)
    r = scan_float(s, 1, n)
    (r === nothing || r[2] != n + 1) && return nothing
    return r[1]
end

# `s` as an Int if it is an optional sign followed by at most 18 digits, else `nothing`.
@inline function _fast_int(s::SubString{String})
    n = ncodeunits(s)
    n == 0 && return nothing
    i = 1
    b = codeunit(s, i)
    neg = b == UInt8('-')
    (neg || b == UInt8('+')) && (i += 1)
    i > n && return nothing
    n - i + 1 > 18 && return nothing
    v = 0
    @inbounds while i <= n
        b = codeunit(s, i)
        isdigit_b(b) || return nothing
        v = v * 10 + Int(b - 0x30)
        i += 1
    end
    return neg ? -v : v
end

# `parse(T, s)`, reading Float64 and Int without `Base.parse` when `s` is in the
# form `_fast_float` or `_fast_int` reads.
@inline fast_parse(::Type{T}, s::SubString{String}) where {T} = parse(T, s)
@inline function fast_parse(::Type{Float64}, s::SubString{String})
    v = _fast_float(s)
    return v === nothing ? parse(Float64, s) : v
end
@inline function fast_parse(::Type{Int}, s::SubString{String})
    v = _fast_int(s)
    return v === nothing ? parse(Int, s) : v
end

# The token at `start` parsed as a `T`, and the position after it, as returned by
# `iter_ws`. `start == 0` means the previous field ended the line, so the row has
# fewer fields than its section declares; this raises the error `parse` gives on an
# empty token.
@noinline function parse_token(::Type{T}, ws::WordedString, start::Int) where {T}
    start < 1 && return (fast_parse(T, NULL_VIEW), 0)
    tok, nxt = iter_ws(ws, start)
    return (fast_parse(T, tok), nxt)
end

@inline iter_num(::Type{T}, ws::WordedString, start::Int) where {T} = parse_token(T, ws, start)

# For Float64, the token is scanned and parsed in one pass. The result is the same
# value and position as `parse_token`: a number `scan_float` reads exactly that ends
# at a delimiter or the end of the line. Anything else goes to `parse_token`.
@inline function iter_num(::Type{Float64}, ws::WordedString, start::Int)
    len = ws.len
    if 1 <= start <= len
        str = ws.s
        i = start
        @inbounds while i <= len && SPACE_B[codeunit(str, i) + 1]
            i += 1
        end
        if i <= len && @inbounds(codeunit(str, i)) != UInt8('%')
            r = scan_float(str, i, len)
            if r !== nothing
                v, j = r
                (j > len || @inbounds(ENDS_B[codeunit(str, j) + 1])) && return (v, j)
            end
        end
    end
    return parse_token(Float64, ws, start)
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

# The line of `s` starting at byte `pos`, and the byte after its newline, or
# `nothing` past the last line. The lines are those of
# `split(replace(s, "\r\n" => "\n"), "\n")`: a `\r` before a newline is dropped,
# and a string ending in a newline has an empty last line.
@inline function next_line(s::String, pos::Int, n::Int)
    if pos > n
        (pos == n + 1 && (n == 0 || codeunit(s, n) == 0x0a)) &&
            return (SubString(s, n + 1, n), n + 2)
        return nothing
    end
    # `findnext` over the code units uses memchr; over the String it decodes characters.
    nl = findnext(==(0x0a), codeunits(s), pos)
    # `thisind`: a line can end in a multibyte character, and a SubString must end
    # at the first byte of its last character.
    nl === nothing && return (SubString(s, pos, thisind(s, n)), n + 1)
    last = nl > pos && codeunit(s, nl - 1) == 0x0d ? nl - 2 : nl - 1
    return (SubString(s, pos, thisind(s, last)), nl + 1)
end

# Byte index of the first non-whitespace character of `line`, or 0 if it is blank.
@inline function first_nonspace(line::SubString{String})
    n = ncodeunits(line)
    i = 1
    @inbounds while i <= n
        b = codeunit(line, i)
        b >= 0x80 && return something(findnext(!isspace, line, i), 0)
        SPACE_B[b + 1] || return i
        i += 1
    end
    return 0
end

# The number of rows in the table whose first line starts at byte `pos`: the lines
# the parsing loop reads as rows, so blank lines and `%` comment lines inside a
# table are not counted.
function get_arr_len(s::String, n::Int, pos::Int)::Int
    count = 0
    p = pos
    while true
        nxt = next_line(s, p, n)
        nxt === nothing && error("Array beginning at byte $pos was not closed")
        line, p = nxt
        k = first_nonspace(line)
        k == 0 && continue
        line[1] == ']' && return count
        codeunit(line, k) == UInt8('%') || (count += 1)
    end
end

# MATPOWER reads a `rateA` of 0 as no limit. Like PowerModels' `calc_thermal_limits!`,
# replace it with a limit the voltage and angle-difference bounds already imply. The
# angle difference is capped at pi, where the voltage difference is largest, so a case
# with MATPOWER's unconstrained -360/360 gets the widest such limit rather than one near 0.
function thermal_limit(r::T, x::T, angmin::T, angmax::T, fr_vmax::T, to_vmax::T) where {T}
    theta_max = min(max(abs(angmin), abs(angmax)), T(pi))
    z2 = r^2 + x^2
    y_mag = iszero(z2) ? zero(T) : inv(sqrt(z2))
    m_vmax = max(fr_vmax, to_vmax)
    c_max = sqrt(fr_vmax^2 + to_vmax^2 - 2 * fr_vmax * to_vmax * cos(theta_max))
    return y_mag * m_vmax * c_max
end

@inbounds @inline @views function parse_matpower_inner(::Type{T}, ::Type{V}, fname :: String, filtered :: Bool) where {T<:Real, V<:AbstractVector}
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
        # `pos` now points at the next line, which is where `get_arr_len` starts.
        line, pos = nxt
        line_len = line.ncodeunits
        line_ind += 1
        k = first_nonspace(line)
        k != 0 && codeunit(line, k) == UInt8('%') && continue
        if in_array && k != 0
            if startswith(line, "];")
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
                f_bus = bus_map[branch_words[1] - bus_offset]
                t_bus = bus_map[branch_words[2] - bus_offset]
                angmin = branch_words[12] / T(180.0) * T(pi)
                angmax = branch_words[13] / T(180.0) * T(pi)
                rate_a = branch_words[6] / baseMVA
                if rate_a <= 0
                    rate_a = thermal_limit(branch_words[3], branch_words[4], angmin, angmax,
                                           bus[f_bus].vmax, bus[t_bus].vmax)
                end
                branch[row_num] = BranchData{T}(
                    row_num,
                    f_bus,
                    t_bus,
                    branch_words[3],
                    branch_words[4],
                    branch_words[5] / T(2.0),
                    branch_words[5] / T(2.0),
                    T(0.0),
                    T(0.0),
                    rate_a,
                    branch_words[7] / baseMVA,
                    branch_words[8] / baseMVA,
                    branch_words[9],
                    (branch_words[10]) / T(180.0) * T(pi),
                    branch_words[11],
                    angmin,
                    angmax,
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

    # `t_idx` was set to `f_idx` plus the declared branch count. If filtering dropped
    # branches, it is reset to `f_idx` plus the remaining count.
    num_declared = num_branch
    num_branch = length(branch)
    fix_t_idx = num_branch != num_declared
    arc = V{ArcData{T}}(undef, num_branch * 2)
    for (i, b) in enumerate(branch)
        if fix_t_idx
            branch[i] = BranchData{T}(
                b.i, b.f_bus, b.t_bus, b.br_r, b.br_x, b.b_fr, b.b_to,
                b.g_fr, b.g_to, b.rate_a, b.rate_b, b.rate_c, b.tap, b.shift,
                b.status, b.angmin, b.angmax, b.f_idx, b.f_idx + num_branch,
            )
        end
        arc[i] = ArcData(i, b.f_bus, b.rate_a)
        arc[i+num_branch] = ArcData(i+num_branch, b.t_bus, b.rate_a)
    end

    return PowerData(version, baseMVA, bus, gen, branch, arc, storage)
end

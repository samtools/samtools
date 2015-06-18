#!/usr/bin/env luajit

function string:split(sep, n)
    local a, start = {}, 1;
    sep = sep or "%s+";
    repeat
        local b, e = self:find(sep, start);
        if b == nil then
            table.insert(a, self:sub(start));
            break
        end
        a[#a+1] = self:sub(start, b - 1);
        start = e + 1;
        if n and #a == n then
            table.insert(a, self:sub(start));
            break
        end
    until start > #self;
    return a;
end

function io.xopen(fn, mode)
    mode = mode or 'r';
    if fn == nil then return io.stdin;
    elseif fn == '-' then return (mode == 'r' and io.stdin) or io.stdout;
    elseif fn:sub(-3) == '.gz' then return (mode == 'r' and io.popen('gzip -dc ' .. fn, 'r')) or io.popen('gzip > ' .. fn, 'w');
    elseif fn:sub(-4) == '.bz2' then return (mode == 'r' and io.popen('bzip2 -dc ' .. fn, 'r')) or io.popen('bgzip2 > ' .. fn, 'w');
    else return io.open(fn, mode) end
end

local eps = {};

function eps.func(fp)
    fp = fp or io.stdout
    fp:write("/C { dup 255 and 255 div exch dup -8 bitshift 255 and 255 div 3 1 roll -16 bitshift 255 and 255 div 3 1 roll setrgbcolor } bind def\n")
    fp:write("/L { 4 2 roll moveto lineto } bind def\n")
    fp:write("/LX { dup 4 -1 roll exch moveto lineto } bind def\n")
    fp:write("/LY { dup 4 -1 roll moveto exch lineto } bind def\n")
    fp:write("/LS { 3 1 roll moveto show } bind def\n")
    fp:write("/RS { dup stringwidth pop 4 -1 roll exch sub 3 -1 roll moveto show } bind def\n")
    fp:write("/B { 4 copy 3 1 roll exch 6 2 roll 8 -2 roll moveto lineto lineto lineto closepath } bind def\n")
end

function eps.font(ft, size, fp)
    fp = fp or io.stdout
    fp:write(string.format('/FS %d def\n', size));
    fp:write('/FS4 FS 4 div def\n');
    fp:write('/' .. ft .. ' findfont FS scalefont setfont\n');
end

local scale = 8;

if #arg == 0 then
    print("Usage: r2plot.lua <in.txt>");
    os.exit(1)
end

local fp = io.xopen(arg[1]);
local n = tonumber(fp:read());

print('%!PS-Adobe-3.0 EPSF-3.0');
print('%%' .. string.format('BoundingBox: -%d -%d %.3f %.3f\n', 10*scale, scale, (n+1)*scale, (n+1)*scale));
print(string.format('%.3f setlinewidth', scale));
print(string.format('/plot { setgray moveto 0 %d rlineto } def', scale));
print(string.format('/plothalf { setgray moveto 0 %.2f rlineto } def', scale/2));
eps.func();
eps.font('Helvetica', scale-1);

local i = 1;
for l in fp:lines() do
    local t = l:split('\t');
    print(string.format("%d %d FS4 add (%s) RS", (i-1)*scale-2, (i-1)*scale, t[1]));
    for j = 2, #t do
        if tonumber(t[j]) > 0.01 then
            print(string.format('%.2f %.2f %.2f plot stroke', (i-1+.5)*scale, (j-2)*scale, 1.-t[j]));
        end
    end
    i = i + 1;
end
for j = 1, 21 do
    print(string.format('%.2f %.2f %.2f plothalf stroke', -8*scale, (j-1) * scale/2, 1.-(j-1)/20));
end
print('showpage');

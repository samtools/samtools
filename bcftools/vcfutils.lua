#!/usr/bin/env luajit

-----------------------------------
-- BEGIN: routines from klib.lua --
-----------------------------------

-- Description: getopt() translated from the BSD getopt(); compatible with the default Unix getopt()
--[[ Example:
	for o, a in os.getopt(arg, 'a:b') do
		print(o, a)
	end
]]--
function os.getopt(args, ostr)
	local arg, place = nil, 0;
	return function ()
		if place == 0 then -- update scanning pointer
			place = 1
			if #args == 0 or args[1]:sub(1, 1) ~= '-' then place = 0; return nil end
			if #args[1] >= 2 then
				place = place + 1
				if args[1]:sub(2, 2) == '-' then -- found "--"
					table.remove(args, 1);
					place = 0
					return nil;
				end
			end
		end
		local optopt = place <= #args[1] and args[1]:sub(place, place) or nil
		place = place + 1;
		local oli = optopt and ostr:find(optopt) or nil
		if optopt == ':' or oli == nil then -- unknown option
			if optopt == '-' then return nil end
			if place > #args[1] then
				table.remove(args, 1);
				place = 0;
			end
			return '?';
		end
		oli = oli + 1;
		if ostr:sub(oli, oli) ~= ':' then -- do not need argument
			arg = nil;
			if place > #args[1] then
				table.remove(args, 1);
				place = 0;
			end
		else -- need an argument
			if place <= #args[1] then  -- no white space
				arg = args[1]:sub(place);
			else
				table.remove(args, 1);
				if #args == 0 then -- an option requiring argument is the last one
					place = 0;
					if ostr:sub(1, 1) == ':' then return ':' end
					return '?';
				else arg = args[1] end
			end
			table.remove(args, 1);
			place = 0;
		end
		return optopt, arg;
	end
end

-- Description: string split
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

-- Description: smart file open
function io.xopen(fn, mode)
	mode = mode or 'r';
	if fn == nil then return io.stdin;
	elseif fn == '-' then return (mode == 'r' and io.stdin) or io.stdout;
	elseif fn:sub(-3) == '.gz' then return (mode == 'r' and io.popen('gzip -dc ' .. fn, 'r')) or io.popen('gzip > ' .. fn, 'w');
	elseif fn:sub(-4) == '.bz2' then return (mode == 'r' and io.popen('bzip2 -dc ' .. fn, 'r')) or io.popen('bgzip2 > ' .. fn, 'w');
	else return io.open(fn, mode) end
end

-- Description: log gamma function
-- Required by: math.lbinom()
-- Reference: AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
function math.lgamma(z)
	local x;
	x = 0.1659470187408462e-06     / (z+7);
	x = x + 0.9934937113930748e-05 / (z+6);
	x = x - 0.1385710331296526     / (z+5);
	x = x + 12.50734324009056      / (z+4);
	x = x - 176.6150291498386      / (z+3);
	x = x + 771.3234287757674      / (z+2);
	x = x - 1259.139216722289      / (z+1);
	x = x + 676.5203681218835      / z;
	x = x + 0.9999999999995183;
	return math.log(x) - 5.58106146679532777 - z + (z-0.5) * math.log(z+6.5);
end

-- Description: regularized incomplete gamma function
-- Dependent on: math.lgamma()
--[[
  Formulas are taken from Wiki, with additional input from Numerical
  Recipes in C (for modified Lentz's algorithm) and AS245
  (http://lib.stat.cmu.edu/apstat/245).
 
  A good online calculator is available at:
 
    http://www.danielsoper.com/statcalc/calc23.aspx
 
  It calculates upper incomplete gamma function, which equals
  math.igamma(s,z,true)*math.exp(math.lgamma(s))
]]--
function math.igamma(s, z, complement)

	local function _kf_gammap(s, z)
		local sum, x = 1, 1;
		for k = 1, 100 do
			x = x * z / (s + k);
			sum = sum + x;
			if x / sum < 1e-14 then break end
		end
		return math.exp(s * math.log(z) - z - math.lgamma(s + 1.) + math.log(sum));
	end

	local function _kf_gammaq(s, z)
		local C, D, f, TINY;
		f = 1. + z - s; C = f; D = 0.; TINY = 1e-290;
		-- Modified Lentz's algorithm for computing continued fraction. See Numerical Recipes in C, 2nd edition, section 5.2
		for j = 1, 100 do
			local d;
			local a, b = j * (s - j), j*2 + 1 + z - s;
			D = b + a * D;
			if D < TINY then D = TINY end
			C = b + a / C;
			if C < TINY then C = TINY end
			D = 1. / D;
			d = C * D;
			f = f * d;
			if math.abs(d - 1) < 1e-14 then break end
		end
		return math.exp(s * math.log(z) - z - math.lgamma(s) - math.log(f));
	end

	if complement then
		return ((z <= 1 or z < s) and 1 - _kf_gammap(s, z)) or _kf_gammaq(s, z);
	else 
		return ((z <= 1 or z < s) and _kf_gammap(s, z)) or (1 - _kf_gammaq(s, z));
	end
end

matrix = {}

-- Description: chi^2 test for contingency tables
-- Dependent on: math.igamma()
function matrix.chi2(a)
	if #a == 2 and #a[1] == 2 then -- 2x2 table
		local x, z
		x = (a[1][1] + a[1][2]) * (a[2][1] + a[2][2]) * (a[1][1] + a[2][1]) * (a[1][2] + a[2][2])
		if x == 0 then return 0, 1, false end
		z = a[1][1] * a[2][2] - a[1][2] * a[2][1]
		z = (a[1][1] + a[1][2] + a[2][1] + a[2][2]) * z * z / x
		return z, math.igamma(.5, .5 * z, true), true
	else -- generic table
		local rs, cs, n, m, N, z = {}, {}, #a, #a[1], 0, 0
		for i = 1, n do rs[i] = 0 end
		for j = 1, m do cs[j] = 0 end
		for i = 1, n do -- compute column sum and row sum
			for j = 1, m do cs[j], rs[i] = cs[j] + a[i][j], rs[i] + a[i][j] end
		end
		for i = 1, n do N = N + rs[i] end
		for i = 1, n do -- compute the chi^2 statistics
			for j = 1, m do
				local E = rs[i] * cs[j] / N;
				z = z + (a[i][j] - E) * (a[i][j] - E) / E
			end
		end
		return z, math.igamma(.5 * (n-1) * (m-1), .5 * z, true), true;
	end
end

---------------------------------
-- END: routines from klib.lua --
---------------------------------


---------------------
-- BEGIN: commands --
---------------------

-- CMD vcf2bgl: convert PL tagged VCF to Beagle input --
function cmd_vcf2bgl()
	if #arg == 0 then
		print("\nUsage: vcf2bgl.lua <in.vcf>")
		print("\nNB: This command finds PL by matching /(\\d+),(\\d+),(\\d+)/.\n");
		os.exit(1)
	end
	
	local lookup = {}
	for i = 0, 10000 do lookup[i] = string.format("%.4f", math.pow(10, -i/10)) end
	
	local fp = io.xopen(arg[1])
	for l in fp:lines() do
		if l:sub(1, 2) == '##' then -- meta lines; do nothing
		elseif l:sub(1, 1) == '#' then -- sample lines
			local t, s = l:split('\t'), {}
			for i = 10, #t do s[#s+1] = t[i]; s[#s+1] = t[i]; s[#s+1] = t[i] end
			print('marker', 'alleleA', 'alleleB', table.concat(s, '\t'))
		else -- data line
			local t = l:split('\t');
			if t[5] ~= '.' and t[5]:find(",") == nil and #t[5] == 1 and #t[4] == 1 then -- biallic SNP
				local x, z = -1, {};
				if t[9]:find('PL') then
					for i = 10, #t do
						local AA, Aa, aa = t[i]:match('(%d+),(%d+),(%d+)')
						AA = tonumber(AA); Aa = tonumber(Aa); aa = tonumber(aa);
						if AA ~= nil then
							z[#z+1] = lookup[AA]; z[#z+1] = lookup[Aa]; z[#z+1] = lookup[aa];
						else z[#z+1] = 1; z[#z+1] = 1; z[#z+1] = 1; end
					end
					print(t[1]..':'..t[2], t[4], t[5], table.concat(z, '\t'))
				elseif t[9]:find('GL') then
					print('Error: not implemented')
					os.exit(1)
				end
			end
		end
	end
	fp:close()
end

-- CMD bgl2vcf: convert Beagle output to VCF
function cmd_bgl2vcf()
	if #arg < 2 then
		print('Usage: bgl2vcf.lua <in.phased> <in.gprobs>')
		os.exit(1)
	end
	
	local fpp = io.xopen(arg[1]);
	local fpg = io.xopen(arg[2]);
	for lg in fpg:lines() do
		local lp = fpp:read()
		local tp, tg, a = lp:split('%s'), lg:split('%s', 4), {}
		if tp[1] == 'I' then
			for i = 3, #tp, 2 do a[#a+1] = tp[i] end
			print('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', table.concat(a, '\t'))
		else
			local chr, pos = tg[1]:match('(%S+):(%d+)$')
			a = {chr, pos, '.', tg[2], tg[3], 30, '.', '.', 'GT'}
			for i = 3, #tp, 2 do
				a[#a+1] = ((tp[i] == tg[2] and 0) or 1) .. '|' .. ((tp[i+1] == tg[2] and 0) or 1)
			end
			print(table.concat(a, '\t'))
		end
	end
	fpg:close(); fpp:close();
end

-- CMD freq: count alleles in each population
function cmd_freq()
	-- parse the command line
	local site_only = true; -- print site allele frequency or not
	for c in os.getopt(arg, 's') do
		if c == 's' then site_only = false end
	end
	if #arg == 0 then
		print("\nUsage: vcfutils.lua freq [-s] <in.vcf> [samples.txt]\n")
		print("NB: 1) This command only considers biallelic variants.")
		print("    2) Apply '-s' to get the allele frequency spectrum.")
		print("    3) 'samples.txt' is TAB-delimited with each line consisting of sample and population.")
		print("")
		os.exit(1)
	end
	
	-- read the sample-population pairs
	local pop, sample = {}, {}
	if #arg > 1 then
		local fp = io.xopen(arg[2]);
		for l in fp:lines() do
			local s, p = l:match("^(%S+)%s+(%S+)"); -- sample, population pair
			sample[s] = p; -- FIXME: check duplications
			if pop[p] then table.insert(pop[p], s)
			else pop[p] = {s} end
		end
		fp:close();
	end
	pop['NA'] = {}
	
	-- parse VCF
	fp = (#arg >= 2 and io.xopen(arg[1])) or io.stdin;
	local col, cnt = {}, {};
	for k in pairs(pop) do
		col[k], cnt[k] = {}, {[0]=0};
	end
	for l in fp:lines() do
		if l:sub(1, 2) == '##' then -- meta lines; do nothing
		elseif l:sub(1, 1) == '#' then -- the sample line
			local t, del_NA = l:split('\t'), true;
			for i = 10, #t do
				local k = sample[t[i]]
				if k == nil then
					k, del_NA = 'NA', false
					table.insert(pop[k], t[i])
				end
				table.insert(col[k], i);
				table.insert(cnt[k], 0);
				table.insert(cnt[k], 0);
			end
			if del_NA then pop['NA'], col['NA'], cnt['NA'] = nil, nil, nil end
		else -- data lines
			local t = l:split('\t');
			if t[5] ~= '.' and t[5]:find(",") == nil then -- biallic
				if site_only == true then io.write(t[1], '\t', t[2], '\t', t[4], '\t', t[5]) end
				for k, v in pairs(col) do
					local ac, an = 0, 0;
					for i = 1, #v do
						local a1, a2 = t[v[i]]:match("^(%d).(%d)");
						if a1 ~= nil then ac, an = ac + a1 + a2, an + 2 end
					end
					if site_only == true then io.write('\t', k, ':', an, ':', ac) end
					if an == #cnt[k] then cnt[k][ac] = cnt[k][ac] + 1 end
				end
				if site_only == true then io.write('\n') end
			end
		end
	end
	fp:close();
	
	-- print
	if site_only == false then
		for k, v in pairs(cnt) do
			io.write(k .. "\t" .. #v);
			for i = 0, #v do io.write("\t" .. v[i]) end
			io.write('\n');
		end
	end
end

function cmd_vcf2chi2()
	if #arg < 3 then
		print("Usage: vcfutils.lua vcf2chi2 <in.vcf> <group1.list> <group2.list>");
		os.exit(1)
	end
	
	local g = {};
	
	-- read the list of groups
	local fp = io.xopen(arg[2]);
	for l in fp:lines() do local x = l:match("^(%S+)"); g[x] = 1 end -- FIXME: check duplicate
	fp:close()
	fp = io.xopen(arg[3]);
	for l in fp:lines() do local x = l:match("^(%S+)"); g[x] = 2 end
	fp:close()
	
	-- process VCF
	fp = io.xopen(arg[1])
	local h = {{}, {}}
	for l in fp:lines() do
		if l:sub(1, 2) == '##' then print(l) -- meta lines; do nothing
		elseif l:sub(1, 1) == '#' then -- sample lines
			local t = l:split('\t');
			for i = 10, #t do
				if g[t[i]] == 1 then table.insert(h[1], i)
				elseif g[t[i]] == 2 then table.insert(h[2], i) end
			end
			while #t > 8 do table.remove(t) end
			print(table.concat(t, "\t"))
		else -- data line
			local t = l:split('\t');
			if t[5] ~= '.' and t[5]:find(",") == nil then -- biallic
				local a = {{0, 0}, {0, 0}}
				for i = 1, 2 do
					for _, k in pairs(h[i]) do
						if t[k]:find("^0.0") then a[i][1] = a[i][1] + 2
						elseif t[k]:find("^1.1") then a[i][2] = a[i][2] + 2
						elseif t[k]:find("^0.1") or t[k]:find("^1.0") then
							a[i][1], a[i][2] = a[i][1] + 1, a[i][2] + 1
						end
					end
				end
				local chi2, p, succ = matrix.chi2(a);
				while #t > 8 do table.remove(t) end
				--print(a[1][1], a[1][2], a[2][1], a[2][2], chi2, p);
				if succ then print(table.concat(t, "\t") .. ";PCHI2=" .. string.format("%.3g", p)
						.. string.format(';AF1=%.4g;AF2=%.4g,%.4g', (a[1][2]+a[2][2]) / (a[1][1]+a[1][2]+a[2][1]+a[2][2]),
						a[1][2]/(a[1][1]+a[1][2]), a[2][2]/(a[2][1]+a[2][2])))
				else print(table.concat(t, "\t")) end
			end
		end
	end
	fp:close()
end

-------------------
-- END: commands --
-------------------


-------------------
-- MAIN FUNCTION --
-------------------

if #arg == 0 then
	print("\nUsage:   vcfutils.lua <command> <arguments>\n")
	print("Command: freq        count biallelic alleles in each population")
	print("         vcf2chi2    compute 1-degree chi-square between two groups of samples")
	print("         vcf2bgl     convert PL annotated VCF to Beagle input")
	print("         bgl2vcf     convert Beagle input to VCF")
	print("")
	os.exit(1)
end

local cmd = arg[1]
table.remove(arg, 1)
if cmd == 'vcf2bgl' then cmd_vcf2bgl()
elseif cmd == 'bgl2vcf' then cmd_bgl2vcf()
elseif cmd == 'freq' then cmd_freq()
elseif cmd == 'vcf2chi2' then cmd_vcf2chi2()
else
	print('ERROR: unknown command "' .. cmd .. '"')
	os.exit(1)
end

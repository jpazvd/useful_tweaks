* Nick Cox Statalist Digest 895
program def tolower 
	version 5.0 
	local varlist "opt ex" 
	parse "`*'" 
	parse "`varlist'", parse(" ") 
	while "`1'" != "" { 
		local newname = lower("`1'") 
		cap: rename `1' `newname' 
		mac shift 
	}
end 

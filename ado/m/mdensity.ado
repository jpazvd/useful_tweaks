*! NJC 2.0.0 5 November 2003
* NJC 1.5.2 21 December 1999 
* NJC 1.5.1 13 October 1999 
program mdensity, rclass sort 
	version 8.0
	syntax varlist [if] [in] [fw aw iw], [ AT(varname) N(integer 50) ///
	BY(varname) MISSing Width(numlist >=0) Range(numlist asc min=2 max=2) ///
	LOG LOGIT a(real 0) b(real 1) Zero BIweight COSine EPanechnikov ///
	GAUssian RECtangle PARzen TRIangle YTItle(str) Intensity ROOT ///
	PLOT(str asis) * ]  
	   
	tokenize `varlist'
	local nvars : word count `varlist'

	if "`by'" != "" & `nvars' > 1 {
		di as err "too many variables specified"
		exit 103
	}

	if "`log'" != "" & "`logit'" != "" {
		di as err "must specify only one of log and logit"
		exit 198
	} 	

	marksample touse, novarlist  
	if "`missing'" == "" {
		if "`by'" != "" markout `touse' `by', strok
	}
	else if "`by'" == "" di as txt "missing only applies with by()" 

	qui count if `touse' 
	if r(N) == 0 error 2000 

	if "`at'" == "" { 
		if "`range'" != "" { 
			local min : word 1 of `range' 
			local max : word 2 of `range' 
		}
		else { 
			su `1' if `touse', meanonly 
			local min = r(min) 
			local max = r(max) 
			
			forval i = 2 / `nvars' { 
				su ``i'' if `touse', meanonly  
				local min = min(`min',r(min)) 
				local max = max(`max',r(max)) 
			}  
		}     
	
		tempvar at
		if `n' <= 1 local n = 50 
		if `n' > _N { 
			local n = _N 
			di as txt "n() set to `n'" 
		}	
		qui range `at' `min' `max' `n' 
		if `nvars' == 1 _crcslbl `at' `1'
		else label var `at' "`varlist'" 
	} 

	if "`log'" != "" {
		tempvar logat 
		qui gen `logat' = log(`at')
		local what "`logat'" 
	}
	else if "`logit'" != "" { 
		tempvar logitat 
	    	qui gen `logitat' = log((`at' -`a')/(`b' - `at')) 
	    	local what `logitat' 
	    	local numer = `b' - `a' 
	} 	
	else local what "`at'" 
	
	if "`width'" == "" local width 0 
	local nw : word count `width' 

	if "`by'" == "" { /* no `by' option */
		if `nw' < `nvars' { 
			if `nw' == 1 local width : di _dup(`nvars') "`width' " 
			else { 
				local nextra = `nvars' - `nw' 
				local extra : di _dup(`nextra') " 0" 
				local width "`width'`extra'" 
			}
		}    

	    	qui forval i = 1/`nvars' {
			tempvar d`i'
			
			if "`log'" != "" { 
				tempvar log`i'  
				gen `log`i'' = log(``i'') 
				local which "`log`i''"
			}
			else if "`logit'" != "" {
				tempvar logit`i' 
				gen `logit`i'' = log((``i'' - `a')/(`b' - ``i'')) 
				local which "`logit`i''" 
			}	
			else local which "``i''" 
			
			local w : word `i' of `width' 
		
			kdensity `which' if `touse' [`weight' `exp'], ///
			`biweight' `cosine' `epanechnikov' `gaussian' ///
			`rectangle' `parzen' `triangle'               ///
			width(`w') nograph gen(`d`i'') at(`what')
	    
			local widths "`widths'`r(width)' " 
			local ns "`ns'`r(n)' " 
			local scales "`scales'`r(scale)' "
			
			if "`log'" != "" replace `d`i'' = `d`i'' / `at' 
			else if "`logit'" != "" { 
				replace `d`i'' = ///
			(`d`i'' * `numer') / ((`at' - `a') * (`b' - `at')) 
			}	
			
			if "`intensity'" != "" { 
		    		count if `touse' & `which' < . 
				replace `d`i'' = `d`i'' *  r(N) 
			}
			if "`zero'" == "" replace `d`i'' = . if `d`i'' == 0 
			if "`root'" != "" replace `d`i'' = sqrt(`d`i'') 
			_crcslbl `d`i'' ``i''
			local dlist "`dlist'`d`i'' "
		}
	}

	else { /* by( ) */
		tempvar group 
		qui bysort `touse' `by': gen byte `group' = _n == 1 if `touse'
		qui replace `group' = sum(`group')
		local ng = `group'[_N]
		local bylab : value label `by'
		local vallab : value label `varlist'
	    
		if `nw' < `ng' { 
			if `nw' == 1 local width : di _dup(`ng') "`width' " 
			else { 
				local nextra = `ng' - `nw' 
				local extra : di _dup(`nextra') " 0" 
				local width "`width'`extra'" 
			}
		}   
	    
		if "`log'" != "" { 
			tempvar logvar  
			qui gen `logvar' = log(`varlist') 
			local which "`logvar'"
		}
		else if "`logit'" != "" { 
			tempvar lgtvar 
			qui gen `lgtvar' = ///
				log((`varlist' - `a')/(`b' - `varlist'))
			local which "`lgtvar'" 
		}	
		else local which "`varlist'" 

		qui count if !`touse'
		local j = 1 + r(N)
		qui forval i = 1/`ng' {
			tempvar d`i'     
			local w : word `i' of `width' 
		
			kdensity `which' if `group' == `i' [`weight' `exp'], ///
		        `biweight' `cosine' `epanechnikov' `gaussian' ///
			`rectangle' `parzen' `triangle' ///
			width(`w') nograph gen(`d`i'') at(`what')
	    
			local widths "`widths'`r(width)' " 	
			local ns "`ns'`r(n)' " 
			local scales "`scales'`r(scale)' "
			
			if "`log'" != "" replace `d`i'' = `d`i'' / `at' 
			else if "`logit'" != "" { 
		    		replace `d`i'' = ///
			(`d`i'' * `numer') / ((`at' - `a') * (`b' - `at')) 
			}
			
			count if `group' == `i' 
			local n = r(N) 
			if "`intensity'" != "" replace `d`i'' = `d`i'' * `n' 
			if "`zero'" == "" replace `d`i'' = . if `d`i'' == 0 
			if "`root'" != "" replace `d`i'' = sqrt(`d`i'') 
			local dlist "`dlist' `d`i''"
			local byval = `by'[`j']
			if "`bylab'" != "" local byval : label `bylab' `byval'
			label var `d`i'' `"`byval'"'
			if "`vallab'" != "" label val `d`i'' `vallab'
			local j = `j' + `n'
		}
	}

    	if `"`ytitle'"' == `""' {
		local ytitle = ///
		cond("`intensity'" != "", "Intensity", "Density") 
	}

	twoway line `dlist' `at' , yti(`"`ytitle'"') `options' || /// 
	`plot' 
	
	local widths = trim("`widths'")
	local ns = trim("`ns'")
	local scales = trim("`scales'")
	return local widths "`widths'" 	
	return local ns "`ns'" 
	return local scales "`scales'" 
end
   

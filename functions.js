/*  Function to Round Numbers */
function round(num, dec) {
  var factor = Math.pow(10,dec);
  return Math.round(num * factor) / factor;
}

/* Compute variance of d-type effect size */
function vd(d,grp1n,grp2n) {
    return (grp1n+grp2n)/(grp1n*grp2n)+Math.pow(d,2)/(2*(grp1n+grp2n));
}

/* small sample size bias adjustment factor */
function j(n1,n2) {
    return (1 - 3/(4*(n1+n2-2)-1))
}
/* small sample size bias adjustment factor */
function jN(n) {
    return (1 - 3/(4*(n-2)-1))
}

/* Hedges g from Cohen's d */
function gfromd(d,n1,n2) {
    return (d * j(n1,n2))
}

/* Hedges g from Cohen's d */
function gfromdN(d,n) {
    return (d * jN(n))
}

/* Variance for Hedges g from Cohen's d */
function vgfromvd(v,n1,n2) {
    return (v * Math.pow(j(n1,n2),2))
}

/* Variance for Hedges g from Cohen's d */
function vgfromvdN(v,n) {
    return (v * Math.pow(jN(n),2))
}

/* 95% confidence interavl */
function lower_d(d,v) {
  return d - 1.959964*Math.sqrt(v) ;
}
function upper_d(d,v) {
  return d + 1.959964*Math.sqrt(v) ;
}

/* SD pooled */
function sdpooled(grp1sd,grp1n,grp2sd,grp2n) {
    return Math.sqrt((grp1sd*grp1sd*(grp1n-1)+grp2sd*grp2sd*(grp2n-1))/(grp1n+grp2n-2)) ;
}

/* SD from SD based on gains */
function sd_raw(sdg,r) {
   return sdg/Math.sqrt(2*(1-r)) ;
}

/* Fisher's Zr transformation */
function Zr(r) {
    return .5*Math.log((1+r)/(1-r)) ;
}

/* Inverse Fisher's Zr transformation */
function InvZr(z) {
     return (Math.exp(2*z) - 1) / (Math.exp(2*z) + 1) ;
}

/* logged odds ratio */
function loggedOR(a,b,c,d) {
    if (a == 0 || b == 0 || c == 0 || d == 0 ) {
	a = a + .5;
	b = b + .5;
	c = c + .5;
	d = d + .5;
    }
    return Math.log((a*d)/(b*c)) ;
}    
function vloggedOR(a,b,c,d) {
    if (a == 0 || b == 0 || c == 0 || d == 0 ) {
	a = a + .5;
	b = b + .5;
	c = c + .5;
	d = d + .5;
    }
    return 1/a + 1/b + 1/c + 1/d ;
}    

/* Cohen's d; probit method */
function probitD(p1,n1,p2,n2) {
    if (p1 == 0 || p2 == 0) {
        p1 = p1 + .5/n1 ;
	p2 = p2 + .5/n2 ;
    }
    if (p1 == 1 || p2 == 1) {
        p1 = p1 - (.5/n1) ;
	p2 = p2 - (.5/n2) ;
    }
   if (p1<.5)  {
	z1 = ANorm(p1*2)*-1 ;
    }
    else if (p1>.5)  {
	z1 = ANorm((1-p1)/.5) ;
    }
    else  {
	z1=0 ;
    }
    if (p2<.5) {
	z2 = ANorm(p2*2)*-1 ;
    }
    else if (p2>.5) {
	z2 = ANorm((1-p2)/.5) ;
    }
    else {
	z2=0 ;
    } 
    return z1 - z2 ;
}

function vprobitD(p1,n1,p2,n2) {
    if (p1 == 0 || p2 == 0) {
        p1 = p1 + .5/n1 ;
	p2 = p2 + .5/n2 ;
    }
    if (p1 == 1 || p2 == 1) {
        p1 = p1 - .5/n1 ;
	p2 = p2 - .5/n2 ;
    }
    if (p1<.5)  {
	z1 = ANorm(p1*2)*-1 ;
    }
    else if (p1>.5)  {
	z1 = ANorm((1-p1)/.5) ;
    }
    else  {
	z1=0 ;
    }
    if (p2<.5) {
	z2 = ANorm(p2*2)*-1 ;
    }
    else if (p2>.5) {
	z2 = ANorm((1-p2)/.5) ;
    }
    else {
	z2=0 ;
    }
    return (2*Math.PI*p1*(1-p1)* Math.exp(z1*z1)/n1)+(2*Math.PI*p2*(1-p2)*Math.exp(z2*z2)/n2) ;
}

function doutput(d,v,se,g,vg,seg) {
    frm.d.value =  round(d,4);
    frm.v.value  = round(v,4);
    frm.lower.value = round(lower_d(d,v),4) ;
    frm.upper.value = round(upper_d(d,v),4) ;
    frm.se.value = round(se,4);
    frm.g.value =  round(g,4);
    frm.lowerg.value = round(lower_d(g,vg),4) ;
    frm.upperg.value = round(upper_d(g,vg),4) ;
    frm.vg.value  = round(vg,4);
    frm.seg.value = round(seg,4);
}

function doutput2(dl,vl,sel,dc,vc,sec,
		  gl,vgl,segl,gc,vgc,segc) {
    frm.dl.value =  round(dl,4);
    frm.dc.value =  round(dc,4);
    frm.lowerl.value = round(lower_d(dl,vl),4) ;
    frm.upperl.value = round(upper_d(dl,vl),4) ;
    frm.lowerc.value = round(lower_d(dc,vc),4) ;
    frm.upperc.value = round(upper_d(dc,vc),4) ;
    frm.vl.value = round(vl,6);
    frm.vc.value = round(vc,6);
    frm.sel.value = round(sel,6);
    frm.sec.value = round(sec,6);
    frm.gl.value =  round(gl,4);
    frm.gc.value =  round(gc,4);
    frm.lowergl.value = round(lower_d(gl,vgl),4) ;
    frm.uppergl.value = round(upper_d(gl,vgl),4) ;
    frm.lowergc.value = round(lower_d(gc,vgc),4) ;
    frm.uppergc.value = round(upper_d(gc,vgc),4) ;
    frm.vgl.value = round(vgl,6);
    frm.vgc.value = round(vgc,6);
    frm.segl.value = round(segl,6);
    frm.segc.value = round(segc,6);
}

function doutput3(dp,vp,sep,gp,vgp,segp) {
    frm.dp.value =  round(dp,4); 1
    frm.lowerp.value = round(lower_d(dp,vp),4) ;
    frm.upperp.value = round(upper_d(dp,vp),4) ;
    frm.vp.value = round(vp,6);
    frm.sep.value = round(sep,6);
    frm.gp.value =  round(gp,4);
    frm.lowergp.value = round(lower_d(gp,vgp),4) ;
    frm.uppergp.value = round(upper_d(gp,vgp),4) ;
    frm.vgp.value = round(vgp,6);
    frm.segp.value = round(segp,6); 
}


function routput(r,zr,vz,sez) {
    frm.r.value = round(r,4);
    frm.lowerr.value = round(InvZr(lower_d(zr,vz)),4) ;
    frm.upperr.value = round(InvZr(upper_d(zr,vz)),4) ;
    frm.zr.value = round(zr,4) ;
    frm.lowerz.value = round(lower_d(zr,vz),4) ;
    frm.upperz.value = round(upper_d(zr,vz),4) ;
    frm.vz.value = round(vz,4) ;
    frm.sez.value = round(sez,4) ;
}

function oroutput(oddsratio,riskratio,logged_or,logged_rr,v_or,v_rr) {
    frm.oddsratio.value =  round(oddsratio,4);
    frm.riskratio.value =  round(riskratio,4);
    frm.logged_or.value =  round(logged_or,4);
    frm.logged_rr.value =  round(logged_rr,4);
    frm.v_or.value =  round(v_or,4);
    frm.v_rr.value =  round(v_rr,4);
    frm.lower_or.value =  round(Math.exp(lower_d(logged_or,v_or)),4);
    frm.upper_or.value =  round(Math.exp(upper_d(logged_or,v_or)),4);
    frm.lower_rr.value =  round(Math.exp(lower_d(logged_rr,v_rr)),4);
    frm.upper_rr.value =  round(Math.exp(upper_d(logged_rr,v_rr)),4);
    frm.se_or.value =  round(Math.sqrt(v_or),4);
    frm.se_rr.value =  round(Math.sqrt(v_rr),4);
}

/* Compute d based on means and standard deviations */
function dMeansSDs() {
    var grp1m	= parseFloat(frm.grp1m.value) ;
    var grp1sd	= parseFloat(frm.grp1sd.value) ;
    var grp1n	= parseFloat(frm.grp1n.value) ;
    var grp2m	= parseFloat(frm.grp2m.value) ;
    var grp2sd	= parseFloat(frm.grp2sd.value) ;
    var grp2n	= parseFloat(frm.grp2n.value) ;
    sd_pooled = sdpooled(grp1sd,grp1n,grp2sd,grp2n) ;
    d   = (grp1m-grp2m) / sd_pooled ;
    g   = gfromd(d,grp1n,grp2n);
    v   = vd(d,grp1n,grp2n);
    vg  = vgfromvd(v,grp1n,grp2n);
    se  = Math.sqrt(v);
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on means and standard errors */
function dMeansSEs() {
    var grp1m  = parseFloat(frm.grp1m.value) ;
    var grp1se = parseFloat(frm.grp1se.value) ;
    var grp1n  = parseFloat(frm.grp1n.value) ;
    var grp2m  = parseFloat(frm.grp2m.value) ;
    var grp2se = parseFloat(frm.grp2se.value) ;
    var grp2n  = parseFloat(frm.grp2n.value) ;
    grp1sd = grp1se*Math.sqrt(grp1n) ;
    grp2sd = grp2se*Math.sqrt(grp2n) ;
    sd_pooled = sdpooled(grp1sd,grp1n,grp2sd,grp2n) ;
    d = (grp1m-grp2m) / sd_pooled ;
    v = vd(d,grp1n,grp2n); 
    se  = Math.sqrt(v);
    g   = gfromd(d,grp1n,grp2n);
    vg  = vgfromvd(v,grp1n,grp2n);
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on means and full sample standard deviations */
function dMeansFullSD() {
    var grp1m = parseFloat(frm.grp1m.value) ;
    var grp1n = parseFloat(frm.grp1n.value) ;
    var grp2m = parseFloat(frm.grp2m.value) ;
    var grp2n = parseFloat(frm.grp2n.value) ;
    var sdy   = parseFloat(frm.sdy.value) ;
    totaln = grp1n+grp2n ;
    sdpooled = Math.sqrt(((totaln-1)/(totaln-2))*Math.pow(sdy,2) - (Math.pow((grp1m-grp2m),2)*grp1n*grp2n)/(totaln*(totaln-2))) ; 
    d = (grp1m-grp2m) / sdpooled ;
    v = vd(d,grp1n,grp2n); 
    se  = Math.sqrt(v);
    g   = gfromd(d,grp1n,grp2n);
    vg  = vgfromvd(v,grp1n,grp2n);
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
    frm.sdpooled.value = round(sdpooled,4);
}

/* Compute d based on means aggregating subgroups */
function dMeansSubgroups() {
    var grp1m  = new Array() ;
    var grp1sd = new Array() ;
    var grp1n  = new Array() ;
    var grp2m  = new Array() ;
    var grp2sd = new Array() ;
    var grp2n  = new Array() ;
    var sdtype = parseFloat(frm.sdtype.value) ;
    grp1m[0] = parseFloat(frm.grp1s0m.value) ;    grp1m[1] = parseFloat(frm.grp1s1m.value) ;
    grp1m[2] = parseFloat(frm.grp1s2m.value) ;    grp1m[3] = parseFloat(frm.grp1s3m.value) ;
    grp1m[4] = parseFloat(frm.grp1s4m.value) ;    grp1m[5] = parseFloat(frm.grp1s5m.value) ;
    grp1m[6] = parseFloat(frm.grp1s6m.value) ;    grp1m[7] = parseFloat(frm.grp1s7m.value) ;
    grp1m[8] = parseFloat(frm.grp1s8m.value) ;
    grp1sd[0] = parseFloat(frm.grp1s0sd.value) ;    grp1sd[1] = parseFloat(frm.grp1s1sd.value) ;
    grp1sd[2] = parseFloat(frm.grp1s2sd.value) ;    grp1sd[3] = parseFloat(frm.grp1s3sd.value) ;
    grp1sd[4] = parseFloat(frm.grp1s4sd.value) ;    grp1sd[5] = parseFloat(frm.grp1s5sd.value) ;
    grp1sd[6] = parseFloat(frm.grp1s6sd.value) ;    grp1sd[7] = parseFloat(frm.grp1s7sd.value) ;
    grp1sd[8] = parseFloat(frm.grp1s8sd.value) ;
    grp1n[0] = parseFloat(frm.grp1s0n.value) ;    grp1n[1] = parseFloat(frm.grp1s1n.value) ;
    grp1n[2] = parseFloat(frm.grp1s2n.value) ;    grp1n[3] = parseFloat(frm.grp1s3n.value) ;
    grp1n[4] = parseFloat(frm.grp1s4n.value) ;    grp1n[5] = parseFloat(frm.grp1s5n.value) ;
    grp1n[6] = parseFloat(frm.grp1s6n.value) ;    grp1n[7] = parseFloat(frm.grp1s7n.value) ;
    grp1n[8] = parseFloat(frm.grp1s8n.value) ;
    grp2m[0] = parseFloat(frm.grp2s0m.value) ;    grp2m[1] = parseFloat(frm.grp2s1m.value) ;
    grp2m[2] = parseFloat(frm.grp2s2m.value) ;    grp2m[3] = parseFloat(frm.grp2s3m.value) ;
    grp2m[4] = parseFloat(frm.grp2s4m.value) ;    grp2m[5] = parseFloat(frm.grp2s5m.value) ;
    grp2m[6] = parseFloat(frm.grp2s6m.value) ;    grp2m[7] = parseFloat(frm.grp2s7m.value) ;
    grp2m[8] = parseFloat(frm.grp2s8m.value) ;
    grp2sd[0] = parseFloat(frm.grp2s0sd.value) ;    grp2sd[1] = parseFloat(frm.grp2s1sd.value) ;
    grp2sd[2] = parseFloat(frm.grp2s2sd.value) ;    grp2sd[3] = parseFloat(frm.grp2s3sd.value) ;
    grp2sd[4] = parseFloat(frm.grp2s4sd.value) ;    grp2sd[5] = parseFloat(frm.grp2s5sd.value) ;
    grp2sd[6] = parseFloat(frm.grp2s6sd.value) ;    grp2sd[7] = parseFloat(frm.grp2s7sd.value) ;
    grp2sd[8] = parseFloat(frm.grp2s8sd.value) ;
    grp2n[0] = parseFloat(frm.grp2s0n.value) ;    grp2n[1] = parseFloat(frm.grp2s1n.value) ;
    grp2n[2] = parseFloat(frm.grp2s2n.value) ;    grp2n[3] = parseFloat(frm.grp2s3n.value) ;
    grp2n[4] = parseFloat(frm.grp2s4n.value) ;    grp2n[5] = parseFloat(frm.grp2s5n.value) ;
    grp2n[6] = parseFloat(frm.grp2s6n.value) ;    grp2n[7] = parseFloat(frm.grp2s7n.value) ;
    grp2n[8] = parseFloat(frm.grp2s8n.value) ;
    // create some needed variables
    grp1nm = 0 ;    grp2nm = 0 ;    grp1nm2 = 0 ;
    grp2nm2 = 0 ;    grp1totaln = 0 ;    grp2totaln = 0 ;
    grp1v = 0 ;    grp2v = 0 ;    k1 = 0 ;    k2 = 0 ;
    for (i=0;i<=8;i=i+1) {
	if (grp1n[i]>0 && grp1sd[i]>0 && grp1m[i]!="") {
	    k1 = k1 + 1 ;
	    grp1nm = grp1nm + grp1m[i]*grp1n[i] ;
	    grp1nm2 = grp1nm2 + grp1m[i]*grp1m[i]*grp1n[i] ;
	    grp1totaln = grp1totaln + grp1n[i] ;
	    grp1v = grp1v  + (grp1sd[i]*grp1sd[i]*(grp1n[i]-1)) ; 
	}
	if (grp2n[i]>0 && grp2sd[i]>0 && grp2m[i]!="") {
	    k2 = k2 + 1 ;
	    grp2nm = grp2nm + grp2m[i]*grp2n[i] ;
	    grp2nm2 = grp2nm2 + grp2m[i]*grp2m[i]*grp2n[i] ;
	    grp2totaln = grp2totaln + grp2n[i] ;
	    grp2v = grp2v  + (grp2sd[i]*grp2sd[i]*(grp2n[i]-1)) ; 
	}
    } 
    grp1m = grp1nm/grp1totaln ;
    grp2m = grp2nm/grp2totaln ;
    grp1sd =  Math.sqrt(grp1v/(grp1totaln-k1)) ;
    grp2sd =  Math.sqrt(grp2v/(grp2totaln-k2)) ;
    if (sdtype == 2) {
	ssb1 = (grp1nm2 - (grp1nm*grp1nm)/grp1totaln) ;
	ssb2 = (grp2nm2 - (grp2nm*grp2nm)/grp2totaln) ;
	grp1sd =  Math.sqrt((ssb1+grp1v)/(grp1totaln-1)) ;
	grp2sd =  Math.sqrt((ssb2+grp2v)/(grp2totaln-1)) ;
    }
    sd_pooled = Math.sqrt((grp1sd*grp1sd*(grp1totaln-1)+grp2sd*grp2sd*(grp2totaln-1))/(grp1totaln+grp2totaln-2)) ;
    d = (grp1m-grp2m) / sd_pooled ;
    v = vd(d,grp1totaln,grp2totaln) ;
    se  = Math.sqrt(v);
    g   = gfromd(d,grp1totaln,grp2totaln);
    vg  = vd(g,grp1totaln,grp2totaln);
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
    frm.sdpooled.value = round(sd_pooled,4);
}

/* Compute d based on means and standard deviations with pretest */
function dMeansSDsPretest() {
    var sdmethod  = parseFloat(frm.sdmethod.value) ;
    var r         = parseFloat(frm.r.value) ;
    var grp1m	  = parseFloat(frm.grp1m.value) ;
    var grp1sd	  = parseFloat(frm.grp1sd.value) ;
    var grp1n	  = parseFloat(frm.grp1n.value) ;
    var grp2m	  = parseFloat(frm.grp2m.value) ;
    var grp2sd	  = parseFloat(frm.grp2sd.value) ;
    var grp2n	  = parseFloat(frm.grp2n.value) ;
    var grp1mpre  = parseFloat(frm.grp1mpre.value) ;
    var grp1sdpre = parseFloat(frm.grp1sdpre.value) ;
    var grp2mpre  = parseFloat(frm.grp2mpre.value) ;
    var grp2sdpre = parseFloat(frm.grp2sdpre.value) ;
    var t1        = parseFloat(frm.t1.value) ;
    var t2        = parseFloat(frm.t2.value) ;
    if (Math.abs(t1)>0 & Math.abs(t2)>0) {
	grp1r = Math.abs(((grp1m*grp1m*grp1n)-(grp1sdpre*grp1sdpre*t1*t1 + grp1sd*grp1sd*t1*t1))/(2*grp1sdpre*grp1sd*t1*t1)) ;
	grp2r = Math.abs(((grp2m*grp2m*grp2n)-(grp2sdpre*grp2sdpre*t2*t2 + grp2sd*grp2sd*t2*t2))/(2*grp2sdpre*grp2sd*t2*t2)) ;
	r = (grp1r+grp2r)/2 ;
    }
    df = grp1n + grp2n - 2 ;
    if (sdmethod==1) {
	sd_pooled = sdpooled(grp1sd,grp1n,grp2sd,grp2n) ;
	c  =  1 - (3/(4*(grp1n+grp2n-2)-1)) ;
    }
    else if (sdmethod==2) {
	sd_pooled = sdpooled(grp1sdpre,grp1n,grp2sdpre,grp2n) ;
 	c  =  1 - (3/(4*(grp1n+grp2n-2)-1)) ;
    }
    else if (sdmethod==3) {
	sd_pooled = Math.sqrt((grp1sd*grp1sd*(grp1n-1)+grp2sd*grp2sd*(grp2n-1)+grp1sdpre*grp1sdpre*(grp1n-1)+grp2sdpre*grp2sdpre*(grp2n-1))/(2*df)) ;
	c  =  1 - (3/(4*(2*grp1n+2*grp2n-4)-1)) ;
    }
    d   = ((grp1m-grp1mpre)-(grp2m-grp2mpre)) / sd_pooled ;
    g   = c*d ;
    v =  2*(1-r)*((grp1n+grp2n)/(grp1n*grp2n))*((df)/(df-2))*(1+g*g/(2*(1-r)*((grp1n+grp2n)/(grp1n*grp2n))))-(g*g)/(c*c) ;
    vg  =  v*(c*c) ;
    se = Math.sqrt(v) ;
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
}

/* Means from Cluster Randomized Designs (WWC method) */
function dCluster() {
    var grp1m	  = parseFloat(frm.grp1m.value) ;
    var grp2m	  = parseFloat(frm.grp2m.value) ;
    var grp2sd	  = parseFloat(frm.grp2sd.value) ;
    var grp1sd	  = parseFloat(frm.grp1sd.value) ;
    var grp1n	  = parseFloat(frm.grp1n.value) ;
    var grp2n	  = parseFloat(frm.grp2n.value) ;
    var p         = parseFloat(frm.p.value) ;
    var nc        = parseFloat(frm.nc.value) ;
    var vmethod   = parseFloat(frm.vmethod.value) ;
    var seb       = parseFloat(frm.seb.value) ;
    totaln = grp1n + grp2n ;
    n = totaln/nc ;
    h = Math.pow((totaln-2)-2*(n-1)*p,2)/((totaln-2)*Math.pow(1-p,2)+ n*(totaln-2*n)*Math.pow(p,2)+2*(totaln-2*n)*p*(1-p)) ;
    c = 1 - 3/(4*h-1) ;
    lambda = 1 - (2*(n-1)*p)/(totaln - 1) ;
    sd_pooled = sdpooled(grp1sd,grp1n,grp2sd,grp2n) ;
    d = ((grp1m - grp2m)/sd_pooled)*Math.sqrt(lambda) ;
    g = (c*(grp1m - grp2m)/sd_pooled)*Math.sqrt(lambda) ;
    if (vmethod==1) {
	v  = (grp1n+grp2n)/(grp1n*grp2n)*(1+(n-1)*p)+Math.pow(d,2)/(2*h) ;
	vg = Math.pow(c,2)*(grp1n+grp2n)/(grp1n*grp2n)*(1+(n-1)*p)+Math.pow(g,2)/(2*h) ;
    }
    else if (vmethod==2) {
	v   = Math.pow(seb/sd_pooled,2)*lambda + Math.pow(d,2)/(2*h) ;
	vg  = Math.pow(c,2)*Math.pow(seb/sd_pooled,2)*lambda + Math.pow(g,2)/(2*h) ;
    }
    se  = Math.sqrt(v) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}

/* HLM continuous */
function dRegHLM() {
    var b	  = parseFloat(frm.b.value) ;
    var seb	  = parseFloat(frm.seb.value) ;
    var sdy	  = parseFloat(frm.sdy.value) ;
    var grp1n	  = parseFloat(frm.grp1n.value) ;
    var grp2n	  = parseFloat(frm.grp2n.value) ;
    var p         = parseFloat(frm.p.value) ;
    var nc        = parseFloat(frm.nc.value) ;
    totaln = grp1n + grp2n ;
    n = totaln/nc ;
    h = Math.pow((totaln-2)-2*(n-1)*p,2)/((totaln-2)*Math.pow(1-p,2)+ n*(totaln-2*n)*Math.pow(p,2)+2*(totaln-2*n)*p*(1-p)) ;
    c = 1 - 3/(4*h-1) ;
    lambda = 1 - (2*(n-1)*p)/(totaln - 1) ;
    sdpooled = Math.sqrt((sdy*sdy*(totaln - 1)-(b*b*(grp1n*grp2n)/(grp1n+grp2n)))/(totaln - 2)) ;
    d = (b/sdpooled)*Math.sqrt(lambda) ;
    g = (c*b/sdpooled)*Math.sqrt(lambda) ;
    v   = Math.pow(seb/sdpooled,2)*lambda + Math.pow(d,2)/(2*h) ;
    vg  = Math.pow(c,2)*Math.pow(seb/sdpooled,2)*lambda + Math.pow(g,2)/(2*h) ;
    se  = Math.sqrt(v) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ; 
}

/* Compute d based on t-value and unequal sample sizes */
function dtUnequal() {
    var t = parseFloat(frm.t.value) ;
    var grp1n = parseFloat(frm.grp1n.value) ;
    var grp2n = parseFloat(frm.grp2n.value) ;
    d = t*Math.sqrt((grp1n+grp2n)/(grp1n*grp2n));
    v = vd(d,grp1n,grp2n);
    se  = Math.sqrt(v);
    g   = gfromd(d,grp1n,grp2n);
    vg  = vgfromvd(v,grp1n,grp2n);
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on t-value and total sample sizes (assumes samples size
   are equal */
function dtEqual() {
    var t      = parseFloat(frm.t.value) ;
    var totaln = parseFloat(frm.totaln.value) ;
    var grp1n = totaln*.5 ;
    var grp2n = totaln*.5 ;
    d = 2*t/Math.sqrt(totaln);
    v = vd(d,grp1n,grp2n);
    se  = Math.sqrt(v);
    g   = gfromd(d,grp1n,grp2n);
    vg  = vgfromvd(v,grp1n,grp2n);
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on p-value of t and total sample size */
function dtpEqual() {
    var poft = parseFloat(frm.poft.value) ;
    var noft = parseFloat(frm.noft.value) ;
    var grp1n = noft*.5 ;
    var grp2n = noft*.5 ;
    if(noft>5) {
	t = AStudT(poft,noft) ;
	d = 2*t/Math.sqrt(noft);
	v = vd(d,grp1n,grp2n);
	se  = Math.sqrt(v);
	g   = gfromd(d,grp1n,grp2n);
	vg  = vgfromvd(v,grp1n,grp2n);
	seg = Math.sqrt(vg);
    }
    if(noft<=5) {
	d = NaN ; v = NaN ; se = NaN ;
	g = NaN ; vg = NaN ; seg = NaN ;
    }
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on p-value of t and separate sample sizes */
function dtpUnequal() {
    var poft  = parseFloat(frm.poft.value) ;
    var grp1n = parseFloat(frm.grp1n.value) ;
    var grp2n = parseFloat(frm.grp2n.value) ;
    var totaln = grp1n+grp2n ;
    if(totaln>5) {
	t = AStudT(poft,totaln) ;
        d = t*Math.sqrt((grp1n+grp2n)/(grp1n*grp2n));
        v = vd(d,grp1n,grp2n);
	se  = Math.sqrt(v);
	g   = gfromd(d,grp1n,grp2n);
	vg  = vgfromvd(v,grp1n,grp2n);
	seg = Math.sqrt(vg);
    }
    if(totaln<=5) {
	d = NaN ; v = NaN ; se = NaN ;
	g = NaN ; vg = NaN ; seg = NaN ;
    }
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on a one-way F-test (two-groups) and total sample sizes */
function dfEqual() {
    var f      = parseFloat(frm.f.value) ;
    var totaln = parseFloat(frm.totaln.value) ;
    grp1n = totaln*.5 ;
    grp2n = totaln*.5 ;
    d = 2*Math.sqrt((f/totaln));
    v = vd(d,grp1n,grp2n);
    se  = Math.sqrt(v);
    g   = gfromd(d,grp1n,grp2n);
    vg  = vgfromvd(v,grp1n,grp2n);
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on a one-way F-test (two-groups) and unequal sample sizes */
function dfUnequal() {
    var f = parseFloat(frm.f.value) ;
    var grp1n = parseFloat(frm.grp1n.value) ;
    var grp2n = parseFloat(frm.grp2n.value) ;
    d = Math.sqrt((f*(grp1n+grp2n))/(grp1n*grp2n));
    v = vd(d,grp1n,grp2n);
    se  = Math.sqrt(v);
    g   = gfromd(d,grp1n,grp2n);
    vg  = vgfromvd(v,grp1n,grp2n);
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
}

function dOneway3() {
    var frm = document.frm ;
    var fvalue = parseFloat(frm.fvalue.value) ;
    var mean    = new Array() ;
    var n       = new Array() ;
    mean[0] = parseFloat(frm.mean0.value) ;    mean[1] = parseFloat(frm.mean1.value) ;
    mean[2] = parseFloat(frm.mean2.value) ;    mean[3] = parseFloat(frm.mean3.value) ;
    mean[4] = parseFloat(frm.mean4.value) ;    mean[5] = parseFloat(frm.mean5.value) ;
    mean[6] = parseFloat(frm.mean6.value) ;    mean[7] = parseFloat(frm.mean7.value) ;
    n[0] = parseFloat(frm.n0.value) ;    n[1] = parseFloat(frm.n1.value) ;
    n[2] = parseFloat(frm.n2.value) ;    n[3] = parseFloat(frm.n3.value) ;
    n[4] = parseFloat(frm.n4.value) ;    n[5] = parseFloat(frm.n5.value) ;
    n[6] = parseFloat(frm.n6.value) ;    n[7] = parseFloat(frm.n7.value) ;
    var rowtype = new Array() ;
    for (i=0;i<=2;i+=1)  {
	if (frm.rowtype0[i].checked) {rowtype[0]=i ;}
	if (frm.rowtype1[i].checked) {rowtype[1]=i ;}
	if (frm.rowtype2[i].checked) {rowtype[2]=i ;}
	if (frm.rowtype3[i].checked) {rowtype[3]=i ;}
	if (frm.rowtype4[i].checked) {rowtype[4]=i ;}
	if (frm.rowtype5[i].checked) {rowtype[5]=i ;}
	if (frm.rowtype6[i].checked) {rowtype[6]=i ;}
	if (frm.rowtype7[i].checked) {rowtype[7]=i ;} 
    }
    totaln = 0 ;
    meanXn = 0 ;
    mean2Xn = 0 ;
    grp1nm = 0 ;
    grp2nm = 0 ;
    grp1n = 0 ;
    grp2n = 0 ;
    k = 0 ;
    k1 = 0 ;
    k2 = 0 ;
    for (i=0;i<=7;i=i+1) {
	if (mean[i]!="" && n[i]>0) {
	    k       += 1 ;
	    totaln  += n[i] ;
	    meanXn  += mean[i]*n[i] ;
	    mean2Xn += mean[i]*mean[i]*n[i] ;
	    if (rowtype[i]==0) {
		grp1nm += mean[i]*n[i] ;
		grp1n  += n[i] ;
		k1     += 1 ;
	    }
	    else if (rowtype[i]==1) {
		grp2nm += mean[i]*n[i] ;
		grp2n  += n[i] ;
		k2     += 1 ;
	    } 
	}
    } 
    msb =  (mean2Xn - (meanXn*meanXn)/totaln)/ (k-1) ;
    sd_pooled = Math.sqrt(msb/fvalue) ;
    grp1m = grp1nm/grp1n ;
    grp2m = grp2nm/grp2n ;
    d = (grp1m - grp2m)/sd_pooled ;
    v = vd(d,grp1n,grp2n) ;
    se  = Math.sqrt(v);
    g   = gfromd(d,grp1n,grp2n);
    vg  = vgfromvd(v,grp1n,grp2n);
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on ANCOVA */
function dAncova() {
    var frm = document.frm ;
    var grp1m = parseFloat(frm.grp1m.value) ;
    var grp1n = parseFloat(frm.grp1n.value) ;
    var grp2m = parseFloat(frm.grp2m.value) ;
    var grp2n = parseFloat(frm.grp2n.value) ;
    var mserror = parseFloat(frm.mserror.value) ;
    var r       = parseFloat(frm.r.value) ;
    totaln = grp1n + grp2n ;
    dferror = totaln - 2 ;
    sd_pooled = Math.sqrt((mserror/(1-r*r))*((dferror-1)/(dferror-2))) ;
    d = (grp1m-grp2m) / sd_pooled ;
    v = vd(d,grp1n,grp2n) ; 
    se  = Math.sqrt(v);
    g   = gfromd(d,grp1n,grp2n);
    vg  = vgfromvd(v,grp1n,grp2n);
    seg = Math.sqrt(vg);
    frm.d.value = round(d,4);
    frm.lower.value = round(lower_d(d,v),4) ;
    frm.upper.value = round(upper_d(d,v),4) ;
    frm.v.value = round(v,4);
    frm.se.value = round(se,4);
    frm.g.value =  round(g,4);
    frm.lowerg.value = round(lower_d(g,vg),4) ;
    frm.upperg.value = round(upper_d(g,vg),4) ;
    frm.vg.value  = round(vg,4);
    frm.seg.value = round(seg,4);    
}

/* Compute d based on 2-way ANOVA */
function dTwoway() {
    var frm = document.frm ;
    var fa      = parseFloat(frm.fa.value) ;
    var fb      = parseFloat(frm.fb.value) ;
    var fab     = parseFloat(frm.fab.value) ;
    var dfb     = parseFloat(frm.dfb.value) ;
    var mserror = parseFloat(frm.mserror.value) ;
    var grp1m = parseFloat(frm.grp1m.value) ;
    var grp1n = parseFloat(frm.grp1n.value) ;
    var grp2m = parseFloat(frm.grp2m.value) ;
    var grp2n = parseFloat(frm.grp2n.value) ;
    dftotal = (grp1n + grp2n) - 1 ;
    dfa = 1 ;
    dfab = dfa*dfb ;
    dferror = dftotal - (dfa+dfb+dfab) ;
    ssa = fa*mserror*dfa ;
    ssb = fb*mserror*dfb ;
    ssab = fab*mserror*dfab ;
    sserror = mserror*dferror ; 
    // determine which radio button is checked
    if (frm.sdtype[0].checked)  {
	sd_pooled = Math.sqrt((ssab+sserror)/(dfab+dferror)) ;
    }
    else if (frm.sdtype[1].checked)  {
	sd_pooled = Math.sqrt((ssb+ssab+sserror)/(dfb+dfab+dferror)) ;
    }
    d = (grp1m-grp2m) / sd_pooled ;
    v = vd(d,grp1n,grp2n) ; 
     se  = Math.sqrt(v);
    g   = gfromd(d,grp1n,grp2n);
    vg  = vgfromvd(v,grp1n,grp2n);
    seg = Math.sqrt(vg);
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on 2 by 2 frequency data */
function d2by2() {
    var frm = document.frm ;
    var grp1s = parseFloat(frm.grp1s.value) ;
    var grp1f = parseFloat(frm.grp1f.value) ;
    var grp2s = parseFloat(frm.grp2s.value) ;
    var grp2f = parseFloat(frm.grp2f.value) ;
    var grp1n = grp1s+grp1f ;
    var grp2n = grp2s+grp2f ;
    loggedor  =   loggedOR(grp1s,grp1f,grp2s,grp2f) ;
    vloggedor =  vloggedOR(grp1s,grp1f,grp2s,grp2f) ;
    /* logit method */
    dl  =  Math.sqrt(3)/Math.PI * loggedor ;
    vl  =  3/Math.pow(Math.PI,2) * vloggedor ;
    sel = Math.sqrt(vl) ;
    gl  = gfromd(dl,grp1n,grp2n);
    vgl = vgfromvd(vl,grp1n,grp2n) ;
    segl = Math.sqrt(vgl) ;
    /* Cox logit method */
    dc = loggedor/1.65 ;
    vc = vloggedor*.3673095 ;
    sec = Math.sqrt(vc) ;
    gc  = gfromd(dc,grp1n,grp2n);
    vgc = vgfromvd(vc,grp1n,grp2n) ;
    segc = Math.sqrt(vgc) ;
    /* probit method */
    n1 = grp1s + grp1f ;
    n2 = grp2s + grp2f ;
    p1 = grp1s/n1 ;
    p2 = grp2s/n2 ;
    dp = probitD(p1,n1,p2,n2) ; 
    vp = vprobitD(p1,n1,p2,n2) ;
    sep = Math.sqrt(vp) ;
    gp  = gfromd(dp,grp1n,grp2n) ;
    vgp = vgfromvd(vp,grp1n,grp2n) ;
    segp = Math.sqrt(vgp) ;
    /* output */
    doutput2(dl,vl,sel,
	     dc,vc,sec,
	     gl,vgl,segl,
	     gc,vgc,segc) ;
    doutput3(dp,vp,sep,gp,vgp,segp) ;
}

/* Compute d based on 2 by 2 proportion data */
function dProp() {
    var frm = document.frm ;
    var grp1ps = parseFloat(frm.grp1ps.value) ;
    var grp2ps = parseFloat(frm.grp2ps.value) ;
    var grp1n  = parseFloat(frm.grp1n.value) ;
    var grp2n  = parseFloat(frm.grp2n.value) ;
    var a = grp1ps * grp1n ;
    var b = (1-grp1ps) * grp1n ;
    var c = grp2ps * grp2n ;
    var d = (1-grp2ps) * grp2n ;
    loggedor  =   loggedOR(a,b,c,d) ;
    vloggedor =  vloggedOR(a,b,c,d) ;
    /* logit method */
    dl  =  Math.sqrt(3)/Math.PI * loggedor ;
    vl  =  3/Math.pow(Math.PI,2) * vloggedor ;
    sel = Math.sqrt(vl) ;
    gl  = gfromd(dl,grp1n,grp2n);
    vgl = vgfromvd(vl,grp1n,grp2n) ;
    segl = Math.sqrt(vgl) ;
    /* Cox logit method */
    dc = loggedor/1.65 ;
    vc = vloggedor*.3673095 ;
    sec = Math.sqrt(vc) ;
    gc  = gfromd(dc,grp1n,grp2n);
    vgc = vgfromvd(vc,grp1n,grp2n) ;
    segc = Math.sqrt(vgc) ;
    /* probit method */
    n1 = grp1n ;
    n2 = grp2n ;
    p1 = grp1ps ;
    p2 = grp2ps ;
    dp = probitD(p1,n1,p2,n2) ; 
    vp = vprobitD(p1,n1,p2,n2) ;
    sep = Math.sqrt(vp) ;
    gp  = gfromd(dp,grp1n,grp2n) ;
    vgp = vgfromvd(vp,grp1n,grp2n) ;
    segp = Math.sqrt(vgp) ;
    /* output */
    doutput2(dl,vl,sel,
	     dc,vc,sec,
	     gl,vgl,segl,
	     gc,vgc,segc) ;
    doutput3(dp,vp,sep,gp,vgp,segp) ;
}

/* Compute d based on 2 by 2 percent data */
function dPercent() {
    var frm = document.frm ;
    var grp1ps = parseFloat(frm.grp1ps.value) ;
    var grp2ps = parseFloat(frm.grp2ps.value) ;
    var grp1n  = parseFloat(frm.grp1n.value) ;
    var grp2n  = parseFloat(frm.grp2n.value) ;
    var grp1ps = grp1ps/100 ;
    var grp2ps = grp2ps/100 ;
    var a = grp1ps * grp1n ;
    var b = (1-grp1ps) * grp1n ;
    var c = grp2ps * grp2n ;
    var d = (1-grp2ps) * grp2n ;
    loggedor  =   loggedOR(a,b,c,d) ;
    vloggedor =  vloggedOR(a,b,c,d) ;
    /* logit method */
    dl  =  Math.sqrt(3)/Math.PI * loggedor ;
    vl  =  3/Math.pow(Math.PI,2) * vloggedor ;
    sel = Math.sqrt(vl) ;
    gl  = gfromd(dl,grp1n,grp2n);
    vgl = vgfromvd(vl,grp1n,grp2n) ;
    segl = Math.sqrt(vgl) ;
    /* Cox logit method */
    dc = loggedor/1.65 ;
    vc = vloggedor*.3673095 ;
    sec = Math.sqrt(vc) ;
    gc  = gfromd(dc,grp1n,grp2n);
    vgc = vgfromvd(vc,grp1n,grp2n) ;
    segc = Math.sqrt(vgc) ;
    /* probit method */
    n1 = grp1n ;
    n2 = grp2n ;
    p1 = grp1ps ;
    p2 = grp2ps ;
    dp = probitD(p1,n1,p2,n2) ; 
    vp = vprobitD(p1,n1,p2,n2) ;
    sep = Math.sqrt(vp) ;
    gp  = gfromd(dp,grp1n,grp2n) ;
    vgp = vgfromvd(vp,grp1n,grp2n) ;
    segp = Math.sqrt(vgp) ;
    /* output */
    doutput2(dl,vl,sel,
	     dc,vc,sec,
	     gl,vgl,segl,
	     gc,vgc,segc) ;
    doutput3(dp,vp,sep,gp,vgp,segp) ;
}
    
/* Compute d based on 2 by 2 frequency data with pretest */
function d2by2Pretest() {
    var frm = document.frm ;
    var r = parseFloat(frm.r.value) ;
    var grp1s = parseFloat(frm.grp1s.value) ;
    var grp1f = parseFloat(frm.grp1f.value) ;
    var grp2s = parseFloat(frm.grp2s.value) ;
    var grp2f = parseFloat(frm.grp2f.value) ;
    var grp1n = grp1s+grp1f ;
    var grp2n = grp2s+grp2f ;
    var grp1spre = parseFloat(frm.grp1spre.value) ;
    var grp1fpre = parseFloat(frm.grp1fpre.value) ;
    var grp2spre = parseFloat(frm.grp2spre.value) ;
    var grp2fpre = parseFloat(frm.grp2fpre.value) ;
    var grp1npre = grp1spre+grp1fpre ;
    var grp2npre = grp2spre+grp2fpre ;
    /* post test */
    loggedor  =   loggedOR(grp1s,grp1f,grp2s,grp2f) ;
    vloggedor =  vloggedOR(grp1s,grp1f,grp2s,grp2f) ;
    /* pre test */
    loggedorpre   =   loggedOR(grp1spre,grp1fpre,grp2spre,grp2fpre) ;
    vloggedorpre  =  vloggedOR(grp1s,grp1f,grp2s,grp2fpre) ;
    /* diff in diff */
    loggedor = loggedor - loggedorpre ;
    vloggedor = vloggedor + vloggedorpre - 2*(r*vloggedor*vloggedorpre) ;
    /* logit method */
    dl  =  Math.sqrt(3)/Math.PI * loggedor ;
    vl  =  3/Math.pow(Math.PI,2) * vloggedor ;
    sel = Math.sqrt(vl) ;
    gl  = gfromd(dl,grp1n,grp2n);
    vgl = vgfromvd(vl,grp1n,grp2n) ;
    segl = Math.sqrt(vgl) ;
    /* Cox logit method */
    dc = loggedor/1.65 ;
    vc = vloggedor*.3673095 ;
    sec = Math.sqrt(vc) ;
    gc  = gfromd(dc,grp1n,grp2n);
    vgc = vgfromvd(vc,grp1n,grp2n) ;
    segc = Math.sqrt(vgc) ;
    /* output */
    doutput2(dl,vl,sel,
	     dc,vc,sec,
	     gl,vgl,segl,
	     gc,vgc,segc) ;
}

/* Compute d based on 2 by 2 proportion data with pretest data */
function dPropPretest() {
    var frm = document.frm ;
    var r = parseFloat(frm.r.value) ;
    var grp1ps    = parseFloat(frm.grp1ps.value) ;
    var grp2ps    = parseFloat(frm.grp2ps.value) ;
    var grp1pspre = parseFloat(frm.grp1pspre.value) ;
    var grp2pspre = parseFloat(frm.grp2pspre.value) ;
    var grp1n     = parseFloat(frm.grp1n.value) ;
    var grp2n     = parseFloat(frm.grp2n.value) ;
    var grp1s = grp1ps * grp1n ;
    var grp1f = (1-grp1ps) * grp1n ;
    var grp2s = grp2ps * grp2n ;
    var grp2f = (1-grp2ps) * grp2n ;
    var grp1spre = grp1pspre * grp1n ;
    var grp1fpre = (1-grp1pspre) * grp1n ;
    var grp2spre = grp2pspre * grp2n ;
    var grp2fpre = (1-grp2pspre) * grp2n ;
    /* post test */
    loggedor  =   loggedOR(grp1s,grp1f,grp2s,grp2f) ;
    vloggedor =  vloggedOR(grp1s,grp1f,grp2s,grp2f) ;
    /* pre test */
    loggedorpre   =   loggedOR(grp1spre,grp1fpre,grp2spre,grp2fpre) ;
    vloggedorpre  =  vloggedOR(grp1s,grp1f,grp2s,grp2fpre) ;
    /* diff in diff */
    loggedor = loggedor - loggedorpre ;
    vloggedor = vloggedor + vloggedorpre - 2*(r*vloggedor*vloggedorpre) ;
    /* logit method */
    dl  =  Math.sqrt(3)/Math.PI * loggedor ;
    vl  =  3/Math.pow(Math.PI,2) * vloggedor ;
    sel = Math.sqrt(vl) ;
    gl  = gfromd(dl,grp1n,grp2n);
    vgl = vgfromvd(vl,grp1n,grp2n) ;
    segl = Math.sqrt(vgl) ;
    /* Cox logit method */
    dc = loggedor/1.65 ;
    vc = vloggedor*.3673095 ;
    sec = Math.sqrt(vc) ;
    gc  = gfromd(dc,grp1n,grp2n);
    vgc = vgfromvd(vc,grp1n,grp2n) ;
    segc = Math.sqrt(vgc) ;
    /* output */
    doutput2(dl,vl,sel,
	     dc,vc,sec,
	     gl,vgl,segl,
	     gc,vgc,segc) ;
}

/* Compute d based on 2 by 2 proportion data with pretest data */
function dPercentPretest() {
    var r = parseFloat(frm.r.value) ;
    var grp1ps    = parseFloat(frm.grp1ps.value) ;
    var grp2ps    = parseFloat(frm.grp2ps.value) ;
    var grp1pspre = parseFloat(frm.grp1pspre.value) ;
    var grp2pspre = parseFloat(frm.grp2pspre.value) ;
    var grp1n     = parseFloat(frm.grp1n.value) ;
    var grp2n     = parseFloat(frm.grp2n.value) ;
    var grp1ps    = grp1ps/100 ;
    var grp2ps    = grp2ps/100;
    var grp1pspre = grp1pspre/100 ;
    var grp2pspre = grp2pspre/100 ;
    var grp1s = grp1ps * grp1n ;
    var grp1f = (1-grp1ps) * grp1n ;
    var grp2s = grp2ps * grp2n ;
    var grp2f = (1-grp2ps) * grp2n ;
    var grp1spre = grp1pspre * grp1n ;
    var grp1fpre = (1-grp1pspre) * grp1n ;
    var grp2spre = grp2pspre * grp2n ;
    var grp2fpre = (1-grp2pspre) * grp2n ;
    /* post test */
    loggedor  =   loggedOR(grp1s,grp1f,grp2s,grp2f) ;
    vloggedor =  vloggedOR(grp1s,grp1f,grp2s,grp2f) ;
    /* pre test */
    loggedorpre   =   loggedOR(grp1spre,grp1fpre,grp2spre,grp2fpre) ;
    vloggedorpre  =  vloggedOR(grp1s,grp1f,grp2s,grp2fpre) ;
    /* diff in diff */
    loggedor = loggedor - loggedorpre ;
    vloggedor = vloggedor + vloggedorpre - 2*(r*vloggedor*vloggedorpre) ;
    /* logit method */
    dl  =  Math.sqrt(3)/Math.PI * loggedor ;
    vl  =  3/Math.pow(Math.PI,2) * vloggedor ;
    sel = Math.sqrt(vl) ;
    gl  = gfromd(dl,grp1n,grp2n);
    vgl = vgfromvd(vl,grp1n,grp2n) ;
    segl = Math.sqrt(vgl) ;
    /* Cox logit method */
    dc = loggedor/1.65 ;
    vc = vloggedor*.3673095 ;
    sec = Math.sqrt(vc) ;
    gc  = gfromd(dc,grp1n,grp2n);
    vgc = vgfromvd(vc,grp1n,grp2n) ;
    segc = Math.sqrt(vgc) ;
    /* output */
    doutput2(dl,vl,sel,
	     dc,vc,sec,
	     gl,vgl,segl,
	     gc,vgc,segc) ;
}
    
/* Compute d based on point-biserial r, unequal sample sizes */
function dPointBiUnequal() {
    var r      = parseFloat(frm.r.value) ;
    var grp1n  = parseFloat(frm.grp1n.value) ;
    var grp2n  = parseFloat(frm.grp2n.value) ;
    p = grp1n/(grp1n+grp2n) ;
    d = r / Math.sqrt((1-r*r)*(p*(1-p))) ;
    v = vd(d,grp1n,grp2n) ;
    se = Math.sqrt(v) ;
    g = gfromd(d,grp1n,grp2n) ;
    vg = vd(g,grp1n,grp2n) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on point-biserial r, equal sample sizes */
function dPointBiEqual() {
    var r      = parseFloat(frm.r.value) ;
    var totaln = parseFloat(frm.totaln.value) ;
    var grp1n = totaln*.5 ;
    var grp2n = totaln*.5 ;
    d = 2*r / Math.sqrt((1-r*r)) ;
    v = vd(d,grp1n,grp2n) ;
    se = Math.sqrt(v) ;
    g = gfromd(d,grp1n,grp2n) ;
    vg = vd(g,grp1n,grp2n) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on phi coefficient from 2 by 2 */
function dPhi() {
    var r      = parseFloat(frm.r.value) ;
    var totaln = parseFloat(frm.totaln.value) ;
    var grp1n = totaln*.5 ;
    var grp2n = totaln*.5 ;
    d = 2*r / Math.sqrt(1-r*r) ;
    chisq = r*r * totaln ;
    v = (d*d)/chisq ;
    se = Math.sqrt(v) ;
    g = gfromd(d,grp1n,grp2n) ;
    vg = vd(g,grp1n,grp2n) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on p-value from chi-square (2 by 2) or phi  */
function dPhipChip() {
    var pchisq  = parseFloat(frm.pchisq.value) ;
    var totaln = parseFloat(frm.totaln.value) ;
    var grp1n = totaln*.5 ;
    var grp2n = totaln*.5 ;
    chisq = AChiSq(pchisq,1) ;
    d = 2*Math.sqrt(chisq/(totaln-chisq)) ;
    v = (d*d)/chisq ;
    se = Math.sqrt(v) ;
    g = gfromd(d,grp1n,grp2n) ;
    vg = vd(g,grp1n,grp2n) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on chi-square from 2 by 2 */
function dChisq() {
    var chisq  = parseFloat(frm.chisq.value) ;
    var totaln = parseFloat(frm.totaln.value) ;
    var grp1n = totaln*.5 ;
    var grp2n = totaln*.5 ;
    d = 2*Math.sqrt(chisq/(totaln-chisq)) ;
    v = (d*d)/chisq ;
    se = Math.sqrt(v) ;
    g = gfromd(d,grp1n,grp2n) ;
    vg = vd(g,grp1n,grp2n) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d from a frequency distribution */
function dOrdinalFreq() {
    var grp1f = new Array() ;
    var grp2f = new Array() ;
    grp1f[0] = parseFloat(frm.grp1f0.value) ;  grp2f[0] = parseFloat(frm.grp2f0.value) ;
    grp1f[1] = parseFloat(frm.grp1f1.value) ;  grp2f[1] = parseFloat(frm.grp2f1.value) ;
    grp1f[2] = parseFloat(frm.grp1f2.value) ;  grp2f[2] = parseFloat(frm.grp2f2.value) ;
    grp1f[3] = parseFloat(frm.grp1f3.value) ;  grp2f[3] = parseFloat(frm.grp2f3.value) ;
    grp1f[4] = parseFloat(frm.grp1f4.value) ;  grp2f[4] = parseFloat(frm.grp2f4.value) ;
    grp1f[5] = parseFloat(frm.grp1f5.value) ;  grp2f[5] = parseFloat(frm.grp2f5.value) ;
    grp1f[6] = parseFloat(frm.grp1f6.value) ;  grp2f[6] = parseFloat(frm.grp2f6.value) ;
    grp1f[7] = parseFloat(frm.grp1f7.value) ;  grp2f[7] = parseFloat(frm.grp2f7.value) ;
    grp1f[8] = parseFloat(frm.grp1f8.value) ;  grp2f[8] = parseFloat(frm.grp2f8.value) ;
    grp1f[9] = parseFloat(frm.grp1f9.value) ;  grp2f[9] = parseFloat(frm.grp2f9.value) ;
    grp1n = 0 ;  grp2n = 0 ;  grp1fn = 0 ;
    grp2fn = 0 ;  grp1f2n = 0 ;  grp2f2n = 0 ;
    for (i=0;i<=9;i=i+1) {
	if (grp1f[i]>0) { 
	    grp1fn = grp1fn + i*grp1f[i] ;
	    grp1n  = grp1n  + grp1f[i] ;
	    grp1f2n = grp1f2n + i*i*grp1f[i] ;
	}
	if (grp2f[i]>0) {
	    grp2fn = grp2fn + i*grp2f[i] ;
	    grp2n  = grp2n  + grp2f[i] ;
	    grp2f2n = grp2f2n + i*i*grp2f[i] ;
	}
    }
    grp1m = grp1fn/grp1n ;
    grp2m = grp2fn/grp2n ;
    grp1sd = Math.sqrt((grp1f2n - (grp1fn*grp1fn)/grp1n)/(grp1n-1)) ;
    grp2sd = Math.sqrt((grp2f2n - (grp2fn*grp2fn)/grp2n)/(grp2n-1)) ;
    sd_pooled = sdpooled(grp1sd,grp1n,grp2sd,grp2n) ;
    d = (grp1m-grp2m) / sd_pooled ; 
    g = gfromd(d,grp1n,grp2n) ;
    v = vd(d,grp1n,grp2n) ;
    se = Math.sqrt(v) ;
    vg = vd(g,grp1n,grp2n) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d from a frequency distribution */
function dOrdinalProps() {
    var grp1n  = parseFloat(frm.grp1n.value) ;
    var grp2n  = parseFloat(frm.grp2n.value) ;
    var grp1f = new Array() ;
    var grp2f = new Array() ;
    var grp1p = new Array() ;
    var grp2p = new Array() ;
    grp1p[0] = parseFloat(frm.grp1p0.value) ;  grp2p[0] = parseFloat(frm.grp2p0.value) ;
    grp1p[1] = parseFloat(frm.grp1p1.value) ;  grp2p[1] = parseFloat(frm.grp2p1.value) ;
    grp1p[2] = parseFloat(frm.grp1p2.value) ;  grp2p[2] = parseFloat(frm.grp2p2.value) ;
    grp1p[3] = parseFloat(frm.grp1p3.value) ;  grp2p[3] = parseFloat(frm.grp2p3.value) ;
    grp1p[4] = parseFloat(frm.grp1p4.value) ;  grp2p[4] = parseFloat(frm.grp2p4.value) ;
    grp1p[5] = parseFloat(frm.grp1p5.value) ;  grp2p[5] = parseFloat(frm.grp2p5.value) ;
    grp1p[6] = parseFloat(frm.grp1p6.value) ;  grp2p[6] = parseFloat(frm.grp2p6.value) ;
    grp1p[7] = parseFloat(frm.grp1p7.value) ;  grp2p[7] = parseFloat(frm.grp2p7.value) ;
    grp1p[8] = parseFloat(frm.grp1p8.value) ;  grp2p[8] = parseFloat(frm.grp2p8.value) ;
    grp1p[9] = parseFloat(frm.grp1p9.value) ;  grp2p[9] = parseFloat(frm.grp2p9.value) ;
    grp1f[0] = grp1p[0]*grp1n ;  grp2f[0] = grp2p[0]*grp2n ;
    grp1f[1] = grp1p[1]*grp1n ;  grp2f[1] = grp2p[1]*grp2n ;
    grp1f[2] = grp1p[2]*grp1n ;  grp2f[2] = grp2p[2]*grp2n ;
    grp1f[3] = grp1p[3]*grp1n ;  grp2f[3] = grp2p[3]*grp2n ;
    grp1f[4] = grp1p[4]*grp1n ;  grp2f[4] = grp2p[4]*grp2n ;
    grp1f[5] = grp1p[5]*grp1n ;  grp2f[5] = grp2p[5]*grp2n ;
    grp1f[6] = grp1p[6]*grp1n ;  grp2f[6] = grp2p[6]*grp2n ;
    grp1f[7] = grp1p[7]*grp1n ;  grp2f[7] = grp2p[7]*grp2n ;
    grp1f[8] = grp1p[8]*grp1n ;  grp2f[8] = grp2p[8]*grp2n ;
    grp1f[9] = grp1p[9]*grp1n ;  grp2f[9] = grp2p[9]*grp2n ;
    grp1n = 0 ;  grp2n = 0 ;  grp1fn = 0 ;
    grp2fn = 0 ;  grp1f2n = 0 ;  grp2f2n = 0 ;
    for (i=0;i<=9;i=i+1) {
	if (grp1f[i]>0) { 
	    grp1fn = grp1fn + i*grp1f[i] ;
	    grp1n  = grp1n  + grp1f[i] ;
	    grp1f2n = grp1f2n + i*i*grp1f[i] ;
	}
	if (grp2f[i]>0) {
	    grp2fn = grp2fn + i*grp2f[i] ;
	    grp2n  = grp2n  + grp2f[i] ;
	    grp2f2n = grp2f2n + i*i*grp2f[i] ;
	}
    }
    grp1m = grp1fn/grp1n ;
    grp2m = grp2fn/grp2n ;
    grp1sd = Math.sqrt((grp1f2n - (grp1fn*grp1fn)/grp1n)/(grp1n-1)) ;
    grp2sd = Math.sqrt((grp2f2n - (grp2fn*grp2fn)/grp2n)/(grp2n-1)) ;
    sd_pooled = sdpooled(grp1sd,grp1n,grp2sd,grp2n) ;
    d = (grp1m-grp2m) / sd_pooled ;
    g = gfromd(d,grp1n,grp2n) ;
    v = vd(d,grp1n,grp2n) ;
    se = Math.sqrt(v) ;
    vg = vd(g,grp1n,grp2n) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}


/* Compute d from Unstandardized Regression Coefficient */
function dRegUnstand() {
    var vmethod = parseFloat(frm.vmethod.value) ;
    var b      = parseFloat(frm.b.value) ;
    var seb    = parseFloat(frm.seb.value) ;
    var sdy    = parseFloat(frm.sdy.value) ;
    var grp1n  = parseFloat(frm.grp1n.value) ;
    var grp2n  = parseFloat(frm.grp2n.value) ;
    totaln = grp1n + grp2n ;
    t = b/seb ;
    sdpooled = Math.sqrt((sdy*sdy*(totaln - 1)-(b*b*(grp1n*grp2n)/(grp1n+grp2n)))/(totaln - 2)) ;
    d = b/sdpooled ;
    g = gfromd(d,grp1n,grp2n) ;
    if (vmethod == 1) {
	v  = (d/t)*(d/t) ;
	vg = (g/t)*(g/t) ;
    }
    else if (vmethod == 2) {
	v = vd(d,grp1n,grp2n) ;
	vg = vd(g,grp1n,grp2n) 
    }
    se = Math.sqrt(v) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d from Unstandardized Regression Coefficient */
function dRegStand() {
    var vmethod = parseFloat(frm.vmethod.value) ;
    var beta   = parseFloat(frm.beta.value) ;
    var sebeta = parseFloat(frm.sebeta.value) ;
    var grp1n  = parseFloat(frm.grp1n.value) ;
    var grp2n  = parseFloat(frm.grp2n.value) ;
    sdy = 1 ;
    totaln = grp1n+grp2n ;
    t = beta/sebeta ;
    sdx =  Math.sqrt((grp1n - (grp1n*grp1n)/totaln)/totaln) ;
    b = beta/sdx ;
    sdpooled = Math.sqrt((sdy*sdy*(totaln - 1)-(b*b*(grp1n*grp2n)/(grp1n+grp2n)))/(totaln - 2)) ;
    d = b/sdpooled ;
    g = gfromd(d,grp1n,grp2n) ;
    if (vmethod == 1) {
	v  = (d/t)*(d/t) ;
	vg = (g/t)*(g/t) ;
    }
    else if (vmethod == 2) {
	v = vd(d,grp1n,grp2n) ;
	vg = vd(g,grp1n,grp2n) 
    }
    se = Math.sqrt(v) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d from Unstandardized Regression Coefficient */
function dRegLogistic() {
    var loggedor  = parseFloat(frm.b.value) ;
    var sebeta    = parseFloat(frm.seb.value) ;
    var vloggedor = sebeta*sebeta ;
    var grp1n  = parseFloat(frm.grp1n.value) ;
    var grp2n  = parseFloat(frm.grp2n.value) ;
    /* logit method */
    dl  =  Math.sqrt(3)/Math.PI * loggedor ;
    gl  = gfromd(dl,grp1n,grp2n);
    vl  =  3/Math.pow(Math.PI,2) * vloggedor ;
    sel = Math.sqrt(vl) ;
    vgl = vgfromvd(vl,grp1n,grp2n) ;
    segl = Math.sqrt(vgl) ;
    /* Cox logit method */
    dc = loggedor/1.65 ;
    gc  = gfromd(dc,grp1n,grp2n);
    vc = vloggedor*.3673095 ;
    sec = Math.sqrt(vc) ;
    vgc = vgfromvd(vc,grp1n,grp2n) ;
    segc = Math.sqrt(vgc) ;
    /* output */
    doutput2(dl,vl,sel,dc,vc,sec,
		  gl,vgl,segl,gc,vgc,segc)
}

/* Compute d based on gain score means */
function dGainR() {
    var grp1m   = parseFloat(frm.grp1m.value) ;
    var grp2m   = parseFloat(frm.grp2m.value) ;
    var grp1sdg = parseFloat(frm.grp1sd.value) ;
    var grp2sdg = parseFloat(frm.grp2sd.value) ;
    var grp1n   = parseFloat(frm.grp1n.value) ;
    var grp2n   = parseFloat(frm.grp2n.value) ;
    var r1      = parseFloat(frm.r1.value) ;
    var r2      = parseFloat(frm.r2.value) ;
    if (!r1) {
	r1 = .5 ;
    }
    if (!r2) {
	r2 = .5 ;
    }
    r = (r1+r2)/2 ;
    grp1sd = sd_raw(grp1sdg,r1) ;
    grp2sd = sd_raw(grp2sdg,r2) ;
    sd_pooled = sdpooled(grp1sd,grp1n,grp2sd,grp2n) ;
    d = (grp1m-grp2m) / sd_pooled ;
    c  =  1 - (3/(4*(2*grp1n+2*grp2n-4)-1)) ;
    g   = c*d ;
    df = grp1n + grp2n - 2 ;
    v =  2*(1-r)*((grp1n+grp2n)/(grp1n*grp2n))*((df)/(df-2))*(1+g*g/(2*(1-r)*((grp1n+grp2n)/(grp1n*grp2n))))-(g*g)/(c*c) ;
    vg =  v * c*c ;
    se = Math.sqrt(v) ;
    seg = Math.sqrt(vg) ;
    doutput(d,v,se,g,vg,seg) ;
}

/* Compute d based on gain score means */
function dGainT() {
    var sdmethod  = parseFloat(frm.sdmethod.value) ;
    var grp1m = parseFloat(frm.grp1m.value) ;
    var grp2m = parseFloat(frm.grp2m.value) ;
    var grp1n = parseFloat(frm.grp1n.value) ;
    var grp2n = parseFloat(frm.grp2n.value) ;
    var grp1sd1 = parseFloat(frm.grp1sd1.value) ;
    var grp2sd1 = parseFloat(frm.grp2sd1.value) ;
    var grp1sd2 = parseFloat(frm.grp1sd2.value) ;
    var grp2sd2 = parseFloat(frm.grp2sd2.value) ;
    var grp1t = parseFloat(frm.grp1t.value) ;
    var grp2t = parseFloat(frm.grp2t.value) ;
    grp1r = Math.abs(((grp1m*grp1m*grp1n)-(grp1sd1*grp1sd1*grp1t*grp1t + grp1sd2*grp1sd2*grp1t*grp1t))/(2*grp1sd1*grp1sd2*grp1t*grp1t)) ;
    grp2r = Math.abs(((grp2m*grp2m*grp2n)-(grp2sd1*grp2sd1*grp2t*grp2t + grp2sd2*grp2sd2*grp2t*grp2t))/(2*grp2sd1*grp2sd2*grp2t*grp2t)) ;
    if (!grp1r) {
	grp1r = .5 ;
    }
    if (!grp2r) {
	grp2r = .5 ;
    }
    r = (grp1r+grp2r)/2 ;
    df = grp1n + grp2n - 2 ;
    if (sdmethod==1) {
	sd_pooled = sdpooled(grp1sd2,grp1n,grp2sd2,grp2n) ;
	c  =  1 - (3/(4*(grp1n+grp2n-2)-1)) ;
    }
    else if (sdmethod==2) {
	sd_pooled = sdpooled(grp1sd1,grp1n,grp2sd1,grp2n) ;
 	c  =  1 - (3/(4*(grp1n+grp2n-2)-1)) ;
    }
    else if (sdmethod==3) {
	sd_pooled = Math.sqrt((grp1sd2*grp1sd2*(grp1n-1)+grp2sd2*grp2sd2*(grp2n-1)+grp1sd1*grp1sd1*(grp1n-1)+grp2sd1*grp2sd1*(grp2n-1))/(2*df)) ;
	c  =  1 - (3/(4*(2*grp1n+2*grp2n-4)-1)) ;
    } 
    d   = (grp1m-grp2m) / sd_pooled ; 
    g   = c*d ;
    v =  2*(1-r)*((grp1n+grp2n)/(grp1n*grp2n))*((df)/(df-2))*(1+g*g/(2*(1-r)*((grp1n+grp2n)/(grp1n*grp2n))))-(g*g)/(c*c) ;
    vg  =  v*(c*c) ;
    se = Math.sqrt(v) ;
    seg = Math.sqrt(vg); 
    doutput(d,v,se,g,vg,seg) ;
}


/* Compute OR/RR based on 2 by 2 frequency data */
function or(a,b,c,d) {
    if (a == 0 || b == 0 || c == 0 || d == 0 ) {
	a = a + .5;	b = b + .5;
	c = c + .5;	d = d + .5;
    }
    oddsratio = (a*d)/(b*c) ;
    riskratio = (a/(a+b))/(c/(c+d)) ;
    logged_or  = Math.log(oddsratio) ;
    logged_rr  = Math.log(riskratio) ;
    v_or = 1/a + 1/b + 1/c + 1/d ;
    v_rr = ((b/a)/(a+b)) + ((d/c)/(c+d)) ;
    oroutput(oddsratio,riskratio,logged_or,logged_rr,v_or,v_rr) ;
}

/* Compute OR/RR based on 2 by 2 frequency data */
function or2by2() {
    var a = parseFloat(frm.a.value) ;
    var b = parseFloat(frm.b.value) ;
    var c = parseFloat(frm.c.value) ;
    var d = parseFloat(frm.d.value) ;
    or(a,b,c,d) ;
}

/* Compute OR/RR based on 2 by 2 proportion data */
function orProportions() {
    var grp1ps = parseFloat(frm.grp1ps.value) ;
    var grp2ps = parseFloat(frm.grp2ps.value) ;
    var grp1n  = parseFloat(frm.grp1n.value) ;
    var grp2n  = parseFloat(frm.grp2n.value) ;
    a = grp1ps*grp1n ;
    b = (1-grp1ps)*grp1n ;
    c = grp2ps*grp2n ;
    d = (1-grp2ps)*grp2n ;
    or(a,b,c,d) ;
}

/* Compute OR/RR based on 2 by 2 percent data */
function orPercents() {
    var grp1ps = parseFloat(frm.grp1ps.value) ;
    var grp2ps = parseFloat(frm.grp2ps.value) ;
    var grp1n  = parseFloat(frm.grp1n.value) ;
    var grp2n  = parseFloat(frm.grp2n.value) ;
    a = grp1ps*grp1n/100 ;
    b = (100-grp1ps)*grp1n/100 ;
    c = grp2ps*grp2n/100 ;
    d = (100-grp2ps)*grp2n/100 ;
    or(a,b,c,d) ;
}

/* Compute OR/RR based on Cohen's d */
function orCohensD() {
    var d  = parseFloat(frm.d.value) ;
    var se = parseFloat(frm.sed.value) ;
    logged_or1 =  Math.PI/Math.sqrt(3) * d ;
    logged_or2 = d*1.65 ;
    or1 = Math.exp(logged_or1) ;
    or2 = Math.exp(logged_or2) ;
    v_or1 = Math.pow(Math.PI,2)/3 * se*se ;
    v_or2 = (se*se)/.367 ;
    se_or1 = Math.sqrt(v_or1) ;
    se_or2 = Math.sqrt(v_or2) ;
    frm.or1.value =  round(or1,4) ;
    frm.or2.value =  round(or2,4) ;
    frm.lower_or1.value =  round(Math.exp(lower_d(logged_or1,v_or1)),4) ;
    frm.upper_or1.value =  round(Math.exp(upper_d(logged_or1,v_or1)),4) ; 
    frm.lower_or2.value =  round(Math.exp(lower_d(logged_or2,v_or2)),4) ;
    frm.upper_or2.value =  round(Math.exp(upper_d(logged_or2,v_or2)),4) ;
    frm.logged_or1.value =  round(logged_or1,4) ;
    frm.logged_or2.value =  round(logged_or2,4) ;
    frm.v_or1.value =  round(v_or1,4) ;
    frm.v_or2.value =  round(v_or2,4) ;
    frm.se_or1.value =  round(se_or1,4) ;
    frm.se_or2.value =  round(se_or2,4) ;
}

function or2by2Pretest() {
    var r = parseFloat(frm.r.value) ;
    var grp1s = parseFloat(frm.grp1s.value) ;
    var grp1f = parseFloat(frm.grp1f.value) ;
    var grp2s = parseFloat(frm.grp2s.value) ;
    var grp2f = parseFloat(frm.grp2f.value) ;
    var grp1n = grp1s+grp1f ;
    var grp2n = grp2s+grp2f ;
    var grp1spre = parseFloat(frm.grp1spre.value) ;
    var grp1fpre = parseFloat(frm.grp1fpre.value) ;
    var grp2spre = parseFloat(frm.grp2spre.value) ;
    var grp2fpre = parseFloat(frm.grp2fpre.value) ;
    var grp1npre = grp1spre+grp1fpre ;
    var grp2npre = grp2spre+grp2fpre ;
    if (grp1s == 0 || grp1f == 0 || grp2s == 0 || grp2f == 0 ) {
	grp1s = grp1s + .5;
	grp1f = grp1f + .5;
	grp2s = grp2s + .5;
	grp2f = grp2f + .5;
    }
    if (grp1spre == 0 || grp1fpre == 0 || grp2spre == 0 || grp2fpre == 0 ) {
	grp1spre = grp1spre + .5;
	grp1fpre = grp1fpre + .5;
	grp2spre = grp2spre + .5;
	grp2fpre = grp2fpre + .5;
    }
    /* post test */
    oddsratio = (grp1s*grp2f)/(grp2s*grp1f) ;
    loggedor  = Math.log(oddsratio) ;
    vloggedor = 1/grp1s + 1/grp1f + 1/grp2s + 1/grp2f ;
    /* pre test */
    oddsratiopre = (grp1spre*grp2fpre)/(grp2spre*grp1fpre) ;
    loggedorpre  = Math.log(oddsratiopre) ;
    vloggedorpre = 1/grp1spre + 1/grp1fpre + 1/grp2spre + 1/grp2fpre ;
    /* diff in diff */
    loggedor = loggedor - loggedorpre ;
    oddsratio = Math.exp(loggedor) ;
    vloggedor = vloggedor + vloggedorpre - 2*(r*vloggedor*vloggedorpre) ;
    seloggedor = Math.sqrt(vloggedor) ;
    /* output */
    frm.oddsratio.value = round(oddsratio, 4) ;
    frm.logged_or.value = round(loggedor, 4) ;
    frm.lower_or.value  = round(lower_d(loggedor,vloggedor),4) ;
    frm.upper_or.value  = round(upper_d(loggedor,vloggedor),4) ;
    frm.v_or.value      = round(vloggedor, 4) ;
    frm.se_or.value     = round(seloggedor, 4) ;
}

/* Correlation (r) from k by j frequency table */
function rkbyj() {
    var freq = new Array(9) ;
    for (i=0; i<9; i++ ) { 
	freq[i] = new Array(9) ;
    }  
    freq[0][0]  = parseFloat(frm.k0j0.value) ;  freq[0][1]  = parseFloat(frm.k0j1.value) ;
    freq[0][2]  = parseFloat(frm.k0j2.value) ;  freq[0][3]  = parseFloat(frm.k0j3.value) ;
    freq[0][4]  = parseFloat(frm.k0j4.value) ;  freq[0][5]  = parseFloat(frm.k0j5.value) ;
    freq[0][6]  = parseFloat(frm.k0j6.value) ;  freq[0][6]  = parseFloat(frm.k0j7.value) ;
    freq[0][8]  = parseFloat(frm.k0j8.value) ;  freq[1][0]  = parseFloat(frm.k1j0.value) ;
    freq[1][1]  = parseFloat(frm.k1j1.value) ;  freq[1][2]  = parseFloat(frm.k1j2.value) ;
    freq[1][3]  = parseFloat(frm.k1j3.value) ;  freq[1][4]  = parseFloat(frm.k1j4.value) ;
    freq[1][5]  = parseFloat(frm.k1j5.value) ;  freq[1][6]  = parseFloat(frm.k1j6.value) ;
    freq[1][6]  = parseFloat(frm.k1j7.value) ;  freq[1][8]  = parseFloat(frm.k1j8.value) ;
    freq[2][0]  = parseFloat(frm.k2j0.value) ;  freq[2][1]  = parseFloat(frm.k2j1.value) ;
    freq[2][2]  = parseFloat(frm.k2j2.value) ;  freq[2][3]  = parseFloat(frm.k2j3.value) ;
    freq[2][4]  = parseFloat(frm.k2j4.value) ;  freq[2][5]  = parseFloat(frm.k2j5.value) ;
    freq[2][6]  = parseFloat(frm.k2j6.value) ;  freq[2][6]  = parseFloat(frm.k2j7.value) ;
    freq[2][8]  = parseFloat(frm.k2j8.value) ;  freq[3][0]  = parseFloat(frm.k3j0.value) ;
    freq[3][1]  = parseFloat(frm.k3j1.value) ;  freq[3][2]  = parseFloat(frm.k3j2.value) ;
    freq[3][3]  = parseFloat(frm.k3j3.value) ;  freq[3][4]  = parseFloat(frm.k3j4.value) ;
    freq[3][5]  = parseFloat(frm.k3j5.value) ;  freq[3][6]  = parseFloat(frm.k3j6.value) ;
    freq[3][6]  = parseFloat(frm.k3j7.value) ;  freq[3][8]  = parseFloat(frm.k3j8.value) ;
    freq[4][0]  = parseFloat(frm.k4j0.value) ;  freq[4][1]  = parseFloat(frm.k4j1.value) ;
    freq[4][2]  = parseFloat(frm.k4j2.value) ;  freq[4][3]  = parseFloat(frm.k4j3.value) ;
    freq[4][4]  = parseFloat(frm.k4j4.value) ;  freq[4][5]  = parseFloat(frm.k4j5.value) ;
    freq[4][6]  = parseFloat(frm.k4j6.value) ;  freq[4][6]  = parseFloat(frm.k4j7.value) ;
    freq[4][8]  = parseFloat(frm.k4j8.value) ;  freq[5][0]  = parseFloat(frm.k5j0.value) ;
    freq[5][1]  = parseFloat(frm.k5j1.value) ;  freq[5][2]  = parseFloat(frm.k5j2.value) ;
    freq[5][3]  = parseFloat(frm.k5j3.value) ;  freq[5][4]  = parseFloat(frm.k5j4.value) ;
    freq[5][5]  = parseFloat(frm.k5j5.value) ;  freq[5][6]  = parseFloat(frm.k5j6.value) ;
    freq[5][6]  = parseFloat(frm.k5j7.value) ;  freq[5][8]  = parseFloat(frm.k5j8.value) ;
    freq[6][0]  = parseFloat(frm.k6j0.value) ;  freq[6][1]  = parseFloat(frm.k6j1.value) ;
    freq[6][2]  = parseFloat(frm.k6j2.value) ;  freq[6][3]  = parseFloat(frm.k6j3.value) ;
    freq[6][4]  = parseFloat(frm.k6j4.value) ;  freq[6][5]  = parseFloat(frm.k6j5.value) ;
    freq[6][6]  = parseFloat(frm.k6j6.value) ;  freq[6][6]  = parseFloat(frm.k6j7.value) ;
    freq[6][8]  = parseFloat(frm.k6j8.value) ;  freq[7][0]  = parseFloat(frm.k7j0.value) ;
    freq[7][1]  = parseFloat(frm.k7j1.value) ;  freq[7][2]  = parseFloat(frm.k7j2.value) ;
    freq[7][3]  = parseFloat(frm.k7j3.value) ;  freq[7][4]  = parseFloat(frm.k7j4.value) ;
    freq[7][5]  = parseFloat(frm.k7j5.value) ;  freq[7][6]  = parseFloat(frm.k7j6.value) ;
    freq[7][6]  = parseFloat(frm.k7j7.value) ;  freq[7][8]  = parseFloat(frm.k7j8.value) ;
    freq[8][0]  = parseFloat(frm.k8j0.value) ;  freq[8][1]  = parseFloat(frm.k8j1.value) ;
    freq[8][2]  = parseFloat(frm.k8j2.value) ;  freq[8][3]  = parseFloat(frm.k8j3.value) ;
    freq[8][4]  = parseFloat(frm.k8j4.value) ;  freq[8][5]  = parseFloat(frm.k8j5.value) ;
    freq[8][6]  = parseFloat(frm.k8j6.value) ;  freq[8][6]  = parseFloat(frm.k8j7.value) ;
    freq[8][8]  = parseFloat(frm.k8j8.value) ;
    fkj = 0 ; fk = 0 ; fj = 0 ; fk2 = 0 ; fj2 = 0 ; n = 0 ;
    for (k=0; k<=8; k++) {
	for (j=0; j<=8; j++) {
	    if (freq[k][j]>0) {
		fkj += freq[k][j] * k * j ;
		fk  += freq[k][j] * k ;
		fj  += freq[k][j] * j ;
		fk2 += freq[k][j] * k*k ;
		fj2 += freq[k][j] * j*j ; 
		n   += freq[k][j] ;
	    }
	}
    } 
    r = (n*fkj-fk*fj)/Math.sqrt((n*fk2 - fk*fk)*(n*fj2-fj*fj))
    zr = Zr(r) ;
    vz = 1/(n-3) ;
    sez = Math.sqrt(vz) ;
    routput(r,zr,vz,sez) ;
}


/* Correlation (r) from r */
function rfromr() {
    var n   = parseFloat(frm.totaln.value) ; 
    var r   = parseFloat(frm.rin.value) ;
    zr  = Zr(r) ;
    vz  = 1/(n-3) ;
    sez = Math.sqrt(vz) ;
    routput(r,zr,vz,sez) ; 
}

/* Correlation (r) from Means and SDs */
function rMeansSDs() {
    var grp1m = parseFloat(frm.grp1m.value) ;
    var grp1sd = parseFloat(frm.grp1sd.value) ;
    var grp1n = parseFloat(frm.grp1n.value) ;
    var grp2m = parseFloat(frm.grp2m.value) ;
    var grp2sd = parseFloat(frm.grp2sd.value) ;
    var grp2n = parseFloat(frm.grp2n.value) ;
    sd_pooled = sdpooled(grp1sd,grp1n,grp2sd,grp2n) ;
    d = (grp1m-grp2m) / sd_pooled ;
    p = grp1n/(grp1n+grp2n) ;
    r = d / Math.sqrt(d*d + 1/(p*(1-p))) ;
    zr = Zr(r) ; 
    vz  = 1/(grp1n+grp2n-3) ;
    sez = Math.sqrt(vz) ;
    routput(r,zr,vz,sez) ; 
}

/* Correlation (r) from 2 by 2 */
function r2by2() {
    var a = parseFloat(frm.a.value) ;
    var b = parseFloat(frm.b.value) ;
    var c = parseFloat(frm.c.value) ;
    var d = parseFloat(frm.d.value) ;
    r = (a*d - b*c) / Math.sqrt((a+b)*(c+d)*(a+c)*(b+d)) ;
    // base v on a rescaled v from a logged odds-ratio (this is more precise than n-3
    logged_or  = Math.log((a*d)/(b*c)) ;
    v_or = 1/a + 1/b + 1/c + 1/d ;
    zr = Zr(r) ; 
    vz = ((zr*zr)/(logged_or*logged_or)) * v_or
    sez = Math.sqrt(vz) ;
    routput(r,zr,vz,sez) ; 
}

/* chi-square and N */
function rchisquare() {
    var chisq   = parseFloat(frm.chisq.value) ;
    var n       = parseFloat(frm.totaln.value) ; 
    r = Math.sqrt(chisq/n) ;
    zr = Zr(r) ;
    vz  = 1/(n-3) ;
    sez = Math.sqrt(vz) ;
    routput(r,zr,vz,sez) ; 
}

/* chi-square and N */
function rchisquarep() {
    var pchisq   = parseFloat(frm.chisqp.value) ;
    var n       = parseFloat(frm.totaln.value) ; 
    chisq = AChiSq(pchisq,1) ;
    r = Math.sqrt(chisq/n) ;
    zr = Zr(r) ;
    vz  = 1/(n-3) ;
    sez = Math.sqrt(vz) ;
    routput(r,zr,vz,sez) ;  
}


/* t-test and N */
function rttest() {
    var n   = parseFloat(frm.totaln.value) ; 
    var t   = parseFloat(frm.tvalue.value) ;
    r = t / Math.sqrt( (t*t) + n - 2) ;
    zr = Zr(r) ;
    vz = 1/(n-3) ;
    sez = Math.sqrt(vz) ;
    routput(r,zr,vz,sez) ;  
}


/* t-test p-value and N */
function rttestp() {
    var poft = parseFloat(frm.poft.value) ;
    var n = parseFloat(frm.totaln.value) ;
    t = AStudT(poft,n) ;
    r = t / Math.sqrt( (t*t) + n - 2) ;
    zr = Zr(r) ;
    vz = 1/(n-3) ;
    sez = Math.sqrt(vz) ;
    routput(r,zr,vz,sez) ;  
}


/*  following code for distribution functions is from John C. Pezzullo */

var Pi=Math.PI; var PiD2=Pi/2; var PiD4=Pi/4; var Pi2=2*Pi
var e=2.718281828459045235; var e10 = 1.105170918075647625
var Deg=180/Pi

function ChiSq(x,n) {
    if(x>1000 | n>1000) { var q=Norm((Math.power(x/n,1/3)+2/(9*n)-1)/Math.sqrt(2/(9*n)))/2; if (x>n) {return q}{return 1-q} }
    var p=Math.exp(-0.5*x); if((n%2)==1) { p=p*Math.sqrt(2*x/Pi) }
    var k=n; while(k>=2) { p=p*x/k; k=k-2 }
    var t=p; var a=n; while(t>1e-15*p) { a=a+2; t=t*x/a; p=p+t }
    return 1-p
    }
function Norm(z) { var q=z*z
    if(Math.abs(z)>7) {return (1-1/q+3/(q*q))*Math.exp(-q/2)/(Math.abs(z)*Math.sqrt(PiD2))} {return ChiSq(q,1) }
    }
function StudT(t,n) {
    t=Math.abs(t); var w=t/Math.sqrt(n); var th=Math.atan(w)
    if(n==1) { return 1-th/PiD2 }
    var sth=Math.sin(th); var cth=Math.cos(th)
    if((n%2)==1)
        { return 1-(th+sth*cth*StatCom(cth*cth,2,n-3,-1))/PiD2 }
        else
        { return 1-sth*StatCom(cth*cth,1,n-3,-1) }
    }
function FishF(f,n1,n2) {
    var x=n2/(n1*f+n2)
    if((n1%2)==0) { return StatCom(1-x,n2,n1+n2-4,n2-2)*Math.pow(x,n2/2) }
    if((n2%2)==0){ return 1-StatCom(x,n1,n1+n2-4,n1-2)*Math.pow(1-x,n1/2) }
    var th=Math.atan(Math.sqrt(n1*f/n2)); var a=th/PiD2; var sth=Math.sin(th); var cth=Math.cos(th)
    if(n2>1) { a=a+sth*cth*StatCom(cth*cth,2,n2-3,-1)/PiD2 }
    if(n1==1) { return 1-a }
    var c=4*StatCom(sth*sth,n2+1,n1+n2-4,n2-2)*sth*Math.pow(cth,n2)/Pi
    if(n2==1) { return 1-a+c/2 }
    var k=2; while(k<=(n2-1)/2) {c=c*k/(k-.5); k=k+1 }
    return 1-a+c
    }
function StatCom(q,i,j,b) {
    var zz=1; var z=zz; var k=i; while(k<=j) { zz=zz*q*k/(k-b); z=z+zz; k=k+2 }
    return z
    }
function ANorm(p) { var v=0.5; var dv=0.5; var z=0
    while(dv>1e-6) { z=1/v-1; dv=dv/2; if(Norm(z)>p) { v=v-dv } else { v=v+dv } }
    return z
    }
function AChiSq(p,n) { var v=0.5; var dv=0.5; var x=0
    while(dv>1e-10) { x=1/v-1; dv=dv/2; if(ChiSq(x,n)>p) { v=v-dv } else { v=v+dv } }
    return x
    }
function AStudT(p,n) { var v=0.5; var dv=0.5; var t=0
    while(dv>1e-6) { t=1/v-1; dv=dv/2; if(StudT(t,n)>p) { v=v-dv } else { v=v+dv } }
    return t
    }
function AFishF(p,n1,n2) { var v=0.5; var dv=0.5; var f=0
    while(dv>1e-10) { f=1/v-1; dv=dv/2; if(FishF(f,n1,n2)>p) { v=v-dv } else { v=v+dv } }
    return f
    }

/* end code from John C. Pezzullo */

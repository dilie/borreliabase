var root1, root2;
var strainData, areaName, areaList, spcData, spcGrpData, spcGrpList=[], geoData, pubData, bioSource, repName, dbLink;
var color_sp, color_area;

var plsList, color_pls, gapCont=60;

var myUrl, type;
var cdhit, locus, symbol, pfamId, idsP, hsvgP, locus_old, paraDnd, orthSeq, spH=16, wTreeP=170, wThumbnail=100, unitW=6.6, color_paraName='#d9fdd9', showAA=1, duration = 2000;
var famData;

var orfData, uvList, contigData, acc_genome={}, conFuse=[], cidList, cid_class={};

var r=30, gids, hUnit=16, wsvgTree=120, yEdge=35, hORF=5, base_tick=yEdge-19, orfTranslate;

var stdOid, stdGid;

var osRoot, ospSeq, ospc;

var transTitle, transData, groupT, rangeT, wBar=80, hUnitT=20, topEdge=26, guidelineT;

$(document).ready(function(){
    $("#main_tabs").tabs();
    
    $.ajax({
        async:false, dataType:"json", url:"js-css/strainData.json",
        success: function(data) {
            strainData=data.strain;
            areaName=data.area;
            spcData=data.spc;
            spcGrpData=data.spcGrp;
            geoData=data.geo;
            pubData=data.pub;
            bioSource=data.biosource;

        areaList = $.merge(Object.keys(areaName),[6]);
        Object.keys(spcGrpData).forEach(function(d){ spcGrpList.push(d*1) });

            var colArr = [0, 1.2, 1.9, 4.1, 5.1, 6.8, 7.8, 9.5];
            color_sp = d3.scaleOrdinal()
                .domain(spcGrpList)
                .range(d3.range(1,9).map(function(d){
                    return d3.hsl(360/9 * colArr[d-1], 1, 0.45)}));

            color_area = d3.scaleOrdinal()
                .domain(areaList)
                .range(["red","orange","#84b837","blue","#00d9d9","darkgray"]);

            strainPlot();
            drawMap();

            dbLink=data.dbLink;
            repName = data.repName;

            plsList = data.repList;
            var supPls=data.supPls, supPlsList=Object.keys(supPls);
            color_pls = d3.scaleOrdinal()
                .domain(supPlsList)
                .range([6,0,8,7].map(function(d){ return d3.schemeCategory10[d]}))

            transTitle=data.transTitle;
            transData = data.trans;
            groupT = data.group_locus;
            rangeT = data.range
        }
    });

    myUrl = 'http://' + window.location.href.split('/')[2] + ':8983/solr/';

    $(document).keyup(function(e){
        $("#newsPara").hide();
        var myEle = document.activeElement;
        var exact;
        if (myEle.id != "searchLoc" && myEle.id != "searchSym" && myEle.id!="searchPfam"){
            $("#suggestion").hide(); return
        }

        if (myEle.id=="searchPfam"){
            $("#suggestion").hide();
            if (e.keyCode=='13'){ queryPfam(myEle.value) }
            return
        }

        if (e.keyCode=='13'){ $("#suggestion").hide(); exact=1 }

        if (myEle.id=='searchLoc'){
            if (locus && locus.toUpperCase()==myEle.value.toUpperCase()){
                $('#searchLoc').val(locus);
                $('#searchSym').val(symbol);
                $('#searchPfam').val(pfamId);
                return
            }
            queryLoc(myEle.value, exact)
        } else {
            querySym(myEle.value, exact)
        }
    });

    $(document).click(function(e) {
        $("#suggestion, #newsPara").hide();
        var myid = e.target.id;
        if (myid=='toSearchL'){ queryLoc($("#searchLoc").val(), 1) }
        if (myid=='toSearchS'){ querySym($("#searchSym").val(), 1) }
        if (myid=='toSearchP'){ queryPfam($("#searchPfam").val()) }
    });

    var tmpAcc = {};
    $.ajax({
        async:false, dataType:"json", url:"js-css/data.json",
        success: function(data){
            contigData=data.contig;
            cidList = data.cidList;
            uvList = data.uvList;

            // get acc_genome
            $.each(contigData, function(acc,obj){
                var geno = obj[0],
                    grp = obj[1];
                if (grp[1]) {
                    conFuse.push(acc);
                    for (var i=0; i<grp.length; i++){ get_acc_g(geno,grp[i],acc) }
                } else {
                    get_acc_g(geno,grp,acc)
                }
            });

            $.each(tmpAcc, function(gid,obj){
                var gg = {};
                if (acc_genome[gid]){ gg = acc_genome[gid] }
                $.each(obj, function(rid, acc){
                    if (acc.length>1) {
                        acc = acc.sort(function(a,b) { return (contigData[a][3]? contigData[a][3] : 0) - ((contigData[b][3]? contigData[b][3] : 0)) })
                    }
                    gg[rid] = acc
                });
                acc_genome[gid] = gg
            })   
        }
    });

    drawTree();
    repliconTbl();

    $("#topScr").scroll(function(){$("#synteny").scrollLeft($("#topScr").scrollLeft())});
    $("#synteny").scroll(function(){$("#topScr").scrollLeft($("#synteny").scrollLeft())});

// set orf colors by cdhit:
    for (var i=0; i<cidList.length; i++) { var j = (i+17)%17; cid_class[cidList[i]] = 'cid' + (j+1) }

    d3.csv("js-css/orf.csv", function(data) { orfData=data });

    function get_acc_g(g,p,a){
        var obj = {};
        if (tmpAcc[g]){ obj = tmpAcc[g] }
        var accs = [];
        if(obj[p]){ accs = obj[p] }
        accs.push(a);
        obj[p] = accs;
        tmpAcc[g] = obj
    }

    $.ajax({
        async:false, dataType:"json", url:"js-css/fam.json",
        success: function(data) { famData = data }
    });

    $('#treeP').css("width", wTreeP+'px');
	$("#seqP").mousemove(function(e) {
		var xMouse = e.pageX - this.offsetLeft - $('#nameP').width() + $("#seqP").scrollLeft()-190;
		var yMouse = e.pageY - this.offsetTop - 151;

	   	var ymax = $("#seqPara").height()-15;
    	$("#gdHp").css("top", (yMouse>0? (yMouse<ymax? yMouse : ymax) : 0) +"px");
    	$("#gdHp").css("width", $("#seqP svg").width()+"px");

	    if (xMouse >= $("#seqP svg").width()-5 || xMouse < 0) {return}
    	$("#gdVp").css("left", xMouse+"px");
        
	    var px = parseInt(xMouse/unitW+1);
	    guidelineP.attr("x1",px/ratioP).attr("x2",px/ratioP)
	});

    $("#selectNtAa").click(function() {
        $(this).html(showAA? "Show Nt" : "Show AA");
        if (showAA){
            d3.select("#seqNt").transition().duration(duration).style("opacity", 0)
            if (!$('#seqAA').length) { aaSeqTbl() }
            d3.select("#seqAA").transition().duration(duration*1.25).delay(duration/4).style("opacity", 1);
            showAA = 0
        } else {
            d3.select("#seqAA").transition().duration(duration).style("opacity", 0)
            d3.select("#seqNt").transition().duration(duration*1.25).delay(duration/4).style("opacity", 1);
            showAA = 1
        }
    });

    $.ajax({
        async:false, dataType:"json", url:"js-css/ospcSeq.json",
        success: function(data) { ospSeq = data; ospcTree() }
    });
	$('#searchIsolate').bind('keypress', function(e){ var code = e.keyCode || e.which; if(code == 13) {searchI()} });
	$('#searchAllele').bind('keypress', function(e){ var code = e.keyCode || e.which; if(code == 13) {searchA()} });

    $(':radio').change(function(){ downSeq(Number($(':radio:checked').val())) });

    loadRepliconMenu();
    addTransTitle();
    addTransContainer();
    draw_transbar(2);
    drawSigLeg();
	$('#repliconMenu').change(function(){ draw_transbar($('#repliconMenu').val()) });
    $("#svgTrans").mousemove(function(e) {
		var yMouse = e.pageY - this.offsetTop - 147;
        if (yMouse<0 || yMouse > $("#guidelineT_area").height()){ yMouse=-50 }
	    guidelineT.attr("y1",yMouse).attr("y2",yMouse)
	});
});


var innerR=116, circleR=2, hLg=17;
function strainPlot(){
    var angle=346,
        gap=2,
        lLabel = 57,
        wPie1 = lLabel + gap*4,
        wPieGeo = 9,
        wPieSpc = 86,
        outerR = innerR + circleR*2 + gap + wPieGeo + wPie1 + wPieSpc,
        marginL = 40,
        marginT = 5,
        wsvg = outerR*2 + marginL,
        hsvg = outerR*2 + marginT;

    $.ajax({
        async:false, dataType:"text", url:"js-css/tree.dnd",
        success: function(data) {
            var parseDND = parseNewick(data);
            root1 = d3.hierarchy(parseDND, function(d) { return d.children });
            root2 = d3.hierarchy(parseDND, function(d) { return d.children });
        }
    });
    var root = root1;
    gids = root.leaves().map(function(d){return d.data.name * 1});

    var svg = d3.select("#strainTree").append("svg").attr("width",wsvg).attr("height",hsvg),
        chart = svg.append("g").attr("transform","translate("+(marginL+outerR)+","+(outerR+marginT)+")");

    var cluster = d3.cluster().size([angle, innerR]).separation(function(a,b){return 1});
    cluster(root);
    setRadius(root, root.data.length=0, innerR/maxLength(root));

    var input = d3.select("#show-length input").on("change", changed),
        timeout = setTimeout(function() { input.property("checked", false).each(changed) }, 2000);

//tree
    var link = chart.append("g")
        .attr("class", "blackLine")
        .selectAll("path").data(root.links()).enter()
        .append("path")
        .each(function(d) { d.target.linkNode = this; })
        .attr("d", linkFun);

    var linkExtension = chart.append("g")
        .attr("class", "grayLine")
        .selectAll("path").data(root.links().filter(function(d){return !d.target.children})).enter()
        .append("path")
        .each(function(d) { d.target.linkExtensionNode = this; })
        .attr("d", linkExtensionFun);

    var areaDot = chart.append("g")
        .selectAll("circle").data(root.leaves()).enter()
        .append("circle")
       .attr("class", function(d){ return 'geo' + strainData[d.data.name].geoId })
        .attr("transform", function(d) {return "rotate(" + (d.x-90) + ")translate(" + (d.radius+circleR) + ",0)"})
        .attr("r", circleR)

//pie1: geography
    var pieSingle = d3.pie().sort(null).endAngle(angle/180*Math.PI).value(1),
        pieMerge=d3.pie().sort(null).endAngle(angle/180*Math.PI).value(function(d){return d[1]});

    var arc;
    base = innerR + circleR*2 + gap;
    arc=d3.arc().innerRadius(base).outerRadius(base+wPieGeo);

    var spG = mergeCt("spId");

    var arcSpc = chart.append("g").selectAll("g").data(pieMerge(spG)).enter();
    var arcMain = chart.append("g").selectAll("g").data(pieSingle(root.leaves())).enter();

    arcMain.append("path")
        .attr("class", "whiteStroke")
        .attr("d", arc)
        .attr("fill", function(d){
            var ge = strainData[d.data.data.name].geoId;
            return color_area(geoData[ge].areaId)
        })
        .on('mouseover', function(d){geoTip(strainData[d.data.data.name].geoId,0)})
        .on('mouseout', function(){geoTip(false)});

// pie2: strain
    base += wPieGeo;
    arc=d3.arc().innerRadius(base).outerRadius(base+wPie1);

    arcMain.append("path").attr("d", arc)
        .attr("class", "pieStrain whiteStroke")
        .attr("id", function(d){ return 'strain'+d.data.data.name})
        .on('mouseover', function(d){strainTip(d.data.data.name)})
        .on('mouseout', function(){strainTip(false)});

        arcMain.append("text")
        .attr("dy", ".31em")
        .attr("transform", function(d) { return "rotate(" + (d.data.x-90) + ")translate(" + (base+gap*2) + ",0)rotate(" + (d.data.x<180? 0 : 180) + ")" })
        .text(function(d) { return strainData[d.data.data.name].name })
        .style("text-anchor", function(d){ return d.data.x<180? "start" : "end"})
        .style("font-size", 11 + "px")
        .on('mouseover', function(d){strainTip(d.data.data.name)})
        .on('mouseout', function(){strainTip(false)})
        .on("click", function(d){openLink('pubmed', strainData[d.data.data.name].cit)})
        .classed("finger", function(d){return strainData[d.data.data.name].cit});

// pie3: ospC
    base += wPie1;
    var wPie2=24;
    arc=d3.arc().innerRadius(base).outerRadius(base+wPie2);

    arcMain.append("path")
        .attr("d", arc)
        .attr("fill", "none")
        .attr("stroke", function(d){
            var spId = strainData[d.data.data.name].spId;
            return spId==139? "white" : 'none'
        });

    base += wPie2/2;
    arcMain.append("text")
        .attr("class", "midAnchor")
        .attr("dy", ".31em")
        .attr("transform", function(d) {
            var c = arc.centroid(d),
                x = c[0], y = c[1],
                h = Math.sqrt(x*x + y*y),
                midAngle = (d.startAngle + d.endAngle)*90/Math.PI;
            return "translate(" + (x/h*base) +  ',' + (y/h*base) +  ") rotate(" + midAngle + ")"        
        })
        .text(function(d) {
            var obj = strainData[d.data.data.name];
            return obj.spId==139 && obj.ospc? obj.ospc : ''
        });

//pie4: SNP group
    base += wPie2/2;
    arc=d3.arc().innerRadius(base).outerRadius(base+wPie2);

    var abcd = mergeCt('snp');
    var arcSnp = chart.append("g").selectAll("g").data(pieMerge(abcd)).enter();

    arcSnp.append("path")
        .attr("d", arc)
        .attr("fill", "none")
        .attr("stroke", function(d){return d.data[0]=='N'? 'none' : "white"});

    base += wPie2/2;
    arcSnp.append("text")
        .attr("class", "midAnchor")
        .attr("transform", function(d) {
            var c = arc.centroid(d),
                x = c[0], y = c[1],
                h = Math.sqrt(x*x + y*y),
                midAngle = (d.startAngle + d.endAngle)*90/Math.PI;
            return "translate(" + (x/h*base) +  ',' + (y/h*base) +  ") rotate(" + midAngle + ")"        
        })
        .attr("dy", ".31em")
        .text(function(d){return d.data[0]=='N'? '' : d.data[0]});

// pie5: species
    base -= wPie2*1.5;
    arc=d3.arc().innerRadius(base).outerRadius(base+wPieSpc);

    arcSpc.append("path")
        .attr("class", "pieSpc whiteStroke")
        .attr("d", arc)
        .attr("fill", function(d){return color_sp(spcData[d.data[0]][1])});

    base += gap*2;
    arcSpc.append("text")
        .attr("class", "ita")
        .attr("transform", function(d) {
            var c = arc.centroid(d),
                x = c[0], y = c[1],
                h = Math.sqrt(x*x + y*y),
                midAngle = (d.startAngle + d.endAngle)*90/Math.PI;
            var toRotate = d.data[1]>6,
                bb = toRotate? base+60 : base;
            return "translate(" + (x/h*bb) +  ',' + (y/h*bb) +  ") rotate(" + (midAngle + (toRotate? 0 : (midAngle>180? 90 : -90))) + ")" 
        })
        .attr("dy", ".31em")
        .style("text-anchor", function(d){
            var midAngle = (d.startAngle + d.endAngle)*90/Math.PI,
                toRotate = d.data[1]>6;
            return toRotate? "middle" : (midAngle<180? "start" : "end")
        })
        .text(function(d){ return spcData[d.data[0]][0] });

// scale
    var scaleX = 64,
        scaleY = 75;
    drawScale(svg, 0.01, innerR, maxLength(root), scaleX, scaleY, "strainScale");

// mark
    var mark = chart.append("g").attr("class","color999")
        .attr("transform", "translate(-6," + (wPieSpc/2-outerR-12) + ")");

    mark.append("text").text("Species");
    mark.append("text").attr("dy", "28px").text("SNP grp");
    mark.append("text").attr("dy", (wPie1+wPie2)/2+7 + "px").text("OC type");
    mark.append("text").attr("dy", (wPie1+wPieSpc)/2+18 + "px").text("Strain");
    mark.append("text").attr("dy", wPie1/2+wPieSpc+wPieGeo/2+7 + "px").text("Geo");

// legend
    var lgSpc = d3.select("#lgSpc").append("svg").attr("width",306).attr("height",spcGrpList.length*hLg)
        .append("g").attr("transform", "translate(4,2)")
        .selectAll("g").data(spcGrpList).enter()
        .append("g").attr("transform", function(d,i) { return "translate(0," + (i*hLg) + ")" });
    lgSpc.append("circle").attr("class", "opacity65").attr("cx", 3.6).attr("cy",3.6).attr("r", 3.6).attr("fill", color_sp);
    lgSpc.append("text").attr("x", 16).attr("dy", ".7em")
        .text(function(d) { return spcGrpData[d] });

    var lgArea = d3.select("#lgArea").append("svg").attr("width",250).attr("height",3*hLg)
        .append("g").attr("transform", "translate(4,2)")
        .selectAll("g").data(areaList).enter()
        .append("g").attr("transform", function(d,i){return "translate(" + (i<3?0:140)+","+(i<3?i : i-3)*hLg+")"});
    lgArea.append("rect")
        .attr("x", 1).attr("y", 1)
        .attr("width", 6).attr("height",6)
        .attr("fill", color_area);
    lgArea.append("text").attr("x", 16).attr("dy", ".7em")
        .text(function(d) { return areaName[d]? areaName[d] : "unknown" });

    function changed(){
        var yes = this.checked;
        clearTimeout(timeout);
        var t = d3.transition().duration(750);
        link.transition(t).attr("d", yes? linkFun : linkConstant);
        linkExtension.transition(t).attr("d", yes? linkExtensionFun : linkExtensionConstant);
        areaDot.transition(t).attr("transform", function(d) {
            return "rotate(" + (d.x-90) + ")translate(" + ((yes? d.radius : d.y)+circleR) + ",0)"});
        d3.select("#strainScale").transition(t).style("opacity", yes? 1 : 0)
    }

    function strainTip(active){
        if (!active){
            $("#straintip").hide();
            d3.selectAll('#map circle').attr("r",circleR).style("stroke","none");
            d3.selectAll('.pieStrain').classed('yellowBack', false)
            return
        }
        var obj = strainData[active];
        d3.select('#map'+obj.geoId).attr("r",circleR*2.5).style("stroke","black");
        d3.select('#strain'+active).classed('yellowBack', true)

        var xMouse = d3.event.pageX-80,
            yMouse = d3.event.pageY-95,
            bioSc, isIta;
        if (obj.source){ bioSc = bioSource[obj.source] }
        if (bioSc){isIta = testIta(bioSc)}

        $("#tipSo").html(bioSc? ((isIta? '<i>' : '') + bioSc + (isIta? '</i>' : '')) : 'Unknown');
        if (obj.cit){
            $("#tipCt").html(obj.cit.map(function(pid){
                var o = pubData[pid];
                return o[0]+ ' et al. (' + o[3] + ') <i>' + o[2] + '</i>'
            }).join("; &nbsp;"))
        } else { $("#tipCt").html('Un-published') }
        $("#straintip").css("left", xMouse+"px").css("top", yMouse + "px").show();
    }

    function mergeCt(field){
        var arr = [], ab;
        root.leaves().forEach(function(d){
            var a = strainData[d.data.name][field];
            if (!a){ a='N' }
            if (!ab || a!=ab){ ab = a; arr.push([ab,1]) }
            else { arr[arr.length-1][1] +=1 }
        });
        return arr
    }
}
function geoTip(active,ismap){
    if (!active){
        $("#geotip").hide();
        d3.selectAll('.pieStrain').classed('yellowBack', false)
        return
    }
    var xMouse = d3.event.pageX-50,
        yMouse = d3.event.pageY-70;
    $("#geotip").css("left", xMouse+"px").css("top", yMouse + "px").html(geoData[active].name).show();
    if (!ismap){ return }
    var here = gids.filter(function(d){ return strainData[d].geoId==active }).map(function(d){ return '#strain'+d}).join(',');
    d3.selectAll(here).classed('yellowBack', true)
}

var wMap=390, hMap=200, mapAct=d3.select(null), mapSvg, mapG, mapPath, mapZoom;
function drawMap(){
    var projection = d3.geoMercator()
                        .translate([wMap/2, hMap/2])
                        .scale(60)
                        .center([0,27]);
    mapPath = d3.geoPath().projection(projection);
 
    mapZoom = d3.zoom().scaleExtent([1,6]).on("zoom", mapZoomed);

    mapSvg = d3.select("#map").append("svg").attr("width",wMap).attr("height",hMap);

    mapSvg.append("rect").attr("id","mapBack").attr("width",wMap).attr("height",hMap).on("click", mapReset);
    mapG = mapSvg.append("g");

    mapSvg.call(mapZoom);

    d3.json("https://unpkg.com/world-atlas@1/world/110m.json", function(error, world) {
        var countries = topojson.feature(world, world.objects.countries).features,
            geoId = gids.map(function(d){ return strainData[d].geoId}).filter(function(d){ return geoData[d].locate});
        
        mapG.selectAll(".country")
            .data(countries).enter()
            .append("path")
            .attr("class", "country")
            .attr("d", mapPath)
            .on('mouseover', function(d){ d3.select(this).classed("country_sel",true)})
            .on('mouseout', function(d){ d3.select(this).classed("country_sel",false)})
            .on("click", mapClick);

        mapG.selectAll("circle")
            .data(geoId).enter()
            .append("circle")
            .attr("id", function(d){ return 'map'+d})
//            .attr("class", "geodot")
            .attr("fill", function(d){return color_area(geoData[d].areaId)})
            .attr("r",circleR)
            .attr("cx", function(d){
                var site = geoData[d].locate,
                    coords = projection([site[1], site[0]]);
                return coords[0]
            })
            .attr("cy", function(d){
                var site = geoData[d].locate,
                    coords = projection([site[1], site[0]]);
                return coords[1]
            })
            .on('mouseover', geoTip)
            .on('mouseout', function(){geoTip(false)})
    })
}
function mapClick(d) {
    if (mapAct.node() === this) return mapReset();
        mapAct.classed("active", false);
        mapAct = d3.select(this).classed("active", true);

        var bounds = mapPath.bounds(d),
            dx = bounds[1][0] - bounds[0][0],
            dy = bounds[1][1] - bounds[0][1],
            x = (bounds[0][0] + bounds[1][0]) / 2,
            y = (bounds[0][1] + bounds[1][1]) / 2,
            scale = Math.max(1, Math.min(8, 0.9/Math.max(dx/wMap, dy/hMap))),
            translate = [wMap/2-scale*x, hMap/2-scale*y];

        mapSvg.transition()
            .duration(750)
            .call(mapZoom.transform, d3.zoomIdentity.translate(translate[0],translate[1]).scale(scale))
}
function mapReset() {
    mapAct.classed("active", false);
    mapAct = d3.select(null);
    mapSvg.transition().duration(1000).call(mapZoom.transform, d3.zoomIdentity)
}

function mapZoomed() {
    mapG.style("stroke-width", 1/d3.event.transform.k + "px");
    mapG.attr("transform", d3.event.transform);
    d3.selectAll('#map circle').attr("r",circleR/d3.event.transform.k)
}


var wsvgRep, wsvgSyn;
function repliconTbl(){
    var wCell = 16, sp_l = 16, hContigs = 76;
    wsvgRep = wCell*plsList.length;
    $('#synteny, #topScr').css("width", (1027-wsvgTree-wsvgRep)+'px');

// rep title
    var svgT = d3.select('#conTitle').append("svg").attr("width", wsvgRep).attr("height", hContigs);
    
    var repTitle = svgT.append("g")
        .selectAll("g").data(plsList).enter()
        .append("g")
        .attr("transform", function(d,i) { return "translate(" + (wCell*i) + ",0)" });

    repTitle.append("rect")
            .attr("width",wCell).attr("height", hContigs)
            .attr("class","fillLess whiteStroke")
            .attr("fill", function(d){ return color_pls(repName[d][1]) });

    repTitle.append("svg:a").attr("xlink:href", "#")
            .append("text")
	    	.attr("transform", "translate(12," + (hContigs-3) + "),rotate(-90)")
            .attr("class", function(d){ return 'pp'+d })
	    	.attr("id", function(d){ return 'plasmid'+d })
	    	.text(function(d){ return repName[d][0] + (repName[d][2]? (' (' + repName[d][2] + ')') : '') })
            .on("click", accFromRep)

// rep table
    var svg = d3.select('#replicon_tbl').append("svg").attr("width", wsvgRep).attr("height", hUnit*gids.length);

    $.each(gids, function(i,gid){
        svg.append("line").attr("x1", 0).attr("x2", wsvgRep).attr("y1", hUnit*(i+1)).attr("y2", hUnit*(i+1));
        var conDot = svg.append("g").attr("class", 'gg'+gid).attr("transform", "translate(0," + (hUnit*(i+1)) + ")");
		var con = acc_genome[gid];
        $.each(plsList, function(j,rid){
            var acc = con[rid];
            if (!acc){ return }
            var ar = acc[0],
                isFuse = $.inArray(ar, conFuse) != -1;
            var t = conDot.append("text")
                    .attr("x", wCell*(j+0.5))
                    .attr("class", 'pp'+rid + (isFuse? (' '+contigData[ar][1].map(function(i){return 'pp'+i}).join(' ')) : ''))
                    .html(isFuse? "&#9643" : "&#9642")
                    .on("mouseover", function(){showPos(gid, rid)})
				    .on("mouseout", function(){hidePos()})
                    .on("click", function(){return accFromRep(rid)});
            if (!/^\d/.test(ar)){t.classed("finger", true).on("dblclick", function(){openLink('ncbi', acc)})}
        })
    })
}
function showPos(gid, rid) {
    d3.selectAll('#treeRep #s' + gid + ', #plasmid'+rid + ', .gg'+gid+' .pp'+rid).classed("darkorange", true);
    var xMouse = d3.event.pageX - 25,
        yMouse = gids.indexOf(gid) * hUnit;
    $("#replicontip").css("left", xMouse+"px").css("top", yMouse+"px").show().html(repName[rid][0])
}
function hidePos() {
    d3.selectAll('#treeRep text, #conTitle text, #replicon_tbl text').classed("darkorange", false);
    $("#replicontip").hide()
}

function accFromRep(grp){
    d3.selectAll('#conTitle text, #replicon_tbl text').classed("magenta", false);
    d3.selectAll('.pp' + grp).classed("magenta", true);
    stdOid=0;
    
    var accs = [];
    gids.forEach(function(d){
        var a = acc_genome[d][grp];
        if (a){ accs.push(a[0]) }
    });
    drawContig(accs,1)
}
function accFromId(cid, oid, g){
    stdOid = cid + '-' + oid;
    stdGid = g? g : 0;
    if (!g) {$('#main_tabs').tabs({ active:1 })}

    var accs = orfData.filter(function(d){ return d.cid==cid && d.orth==oid})
            .map(function(d){
                var con = contigData[d.con],
                    gid = con[0],
                    rid = con[1][0]? con[1][0] : con[1],
                    acc = acc_genome[gid][rid];
                    return acc[0]
            })
            .sort(function(a,b){return gids.indexOf(contigData[a][0]) - gids.indexOf(contigData[b][0])});
    drawContig(accs)
}

var scaleSyn, currAcc;
function drawContig(accs, noAln){
    var svgSyn;
    if (!$('#allORF')[0]){
        svgSyn = d3.select("#synteny").append("svg").attr("height", hUnit*gids.length+yEdge-4);
        svgSyn.append("g").attr("id", "allORF").attr("transform", "translate(0," + yEdge + ")")

        scaleSyn = d3.scaleLinear().domain([0, 100]).range([0, 100/r]);
        drawAxLg()
    } else {
        svgSyn = d3.select("#synteny svg");
        d3.selectAll('#synteny svg .shade').remove()
    }

    var toHide, toAdd;
    if (currAcc) {
        toHide = currAcc.filter(function(d){ return accs.indexOf(d)==-1 });
        toAdd = accs.filter(function(d){ return currAcc.indexOf(d)==-1 });
    } else {
        toHide = [];
        toAdd = accs
    }
    currAcc = accs;

    if (toHide.length){ hideCon(toHide, noAln) }
    if (toAdd.length){ addCon(toAdd) }

    // length of fig:
    var lMax;
    if (noAln){
        lMax = 0;
        $.each(currAcc, function(i,acc){
            var l = $('.con' + acc)[0].attributes['len'].value * 1;
            if (l>lMax){ lMax = l }
        })
    } else {
        var end_len=[], startMax=0, remainMax=0;
        $.each(currAcc, function(i,acc){
            var gid = contigData[acc][0],
                n = gids.indexOf(gid);
            var myid = '#'+gid+'_'+stdOid;

            var myorf = $(myid)[0],
                end =  myorf.attributes['end'].value * 1,
                len =  myorf.attributes['len'].value * 1;
            end_len.push([gid, end, len, acc]);

            if (end > startMax) { startMax = end }
            var remain = $('.con' + acc)[0].attributes['len'].value - end;
            if (remain > remainMax) { remainMax = remain }
        });
        lMax = startMax + remainMax
    }
    
    svgSyn.attr("width", lMax/r+1);
    $('#topS').css("width", (lMax/r+1)+'px');
    
    //translate
    if (noAln){
        accs.forEach(function(d){d3.selectAll('.con'+d)/*.transition().duration(2000)*/.attr("transform", "translate(0)")})
        d3.select("#tick_mark")/*.transition().duration(2000)*/.attr("transform", "translate(0)");
        $("#synteny, #topscr").animate({scrollLeft: 0},1000)
        orfTranslate = 0;
        return
    }
    
    var stdTrans = stdGid? (orfTranslate? orfTranslate[stdGid] : 0) : startMax/r-180;
    
    orfTranslate = {};
    var shadeX=[], shadeY=[];
    end_len.forEach(function(d,i){
        var gid = d[0],
            end = d[1],
            len = d[2],
            acc = d[3],
            n = gids.indexOf(gid);

        var toTrans = Math.round((startMax-end)/r*10)/10;

        if (stdGid && gid==stdGid){$("#synteny, #topScr").scrollLeft(toTrans - stdTrans + $('#synteny').scrollLeft())}
        if (!stdGid && !i){$("#synteny, #topScr").scrollLeft(stdTrans)}
        
        d3.select('.con' + acc).attr("transform", "translate(" + toTrans + ")");

        if (!i) { d3.select("#tick_mark").attr("transform", "translate(" + toTrans + ")") }

        orfTranslate[gid] = toTrans;

        //shade points
		var base = hUnit * n;
	    if (!shadeX.length) {
            shadeX.push(startMax, (startMax-len));
            shadeY.push(base-7, base-7)
        }
	    shadeX.push(startMax-len);
        shadeY.push(base)
    });
    
    //  draw shade
    d3.selectAll('#synteny svg .shade:not(#shadeC)').remove();

    var lastX = shadeX[shadeX.length-1],
        lastY = shadeY[shadeY.length-1];
    shadeX.push(lastX, shadeX[0]);
    shadeY.push(lastY+7, lastY+7);
    var points = $.map(shadeX, function(vx,i){ return (vx)/r + ',' + shadeY[i]}).join(' ');

    if (!$('#shadeC')[0]){
        d3.select('#allORF').append("polygon").attr("id", "shadeC").attr("class","shade").style("stroke-opacity",0)
    }
    d3.select('#shadeC').attr("points", points)
        .style("stroke-opacity", 0)
        .transition().duration(2000).delay(toHide.length? 1000 : 0)
            .style("stroke-opacity", 1);

    d3.selectAll('#conTitle text').classed("magenta", false);

    function drawAxLg(){
        var tick_mark = svgSyn.append("g").attr("id", "tick_mark");
        drawAxis(tick_mark, r, contigData[acc_genome[100][1][0]][2], 200, 5, 1000)

        $('#repInstr').html('<span class="step">ORF: </span><b>Click</b><span class="step"> to align orthologs; </span><b>Double Click</b><span class="step"> for annotation</span>').css("left",677+'px');
        drawLg()
    }
}
function hideCon(cons, noAln){
    cons.forEach(function(d){
        d3.selectAll('.con'+d)
            .transition().duration(noAln? 0 : 2000)
                .attr("transform", function(){return "translate(" + (-$(this).attr("len")/r-1) + "," + 0 + ")"});
        var arr = contigData[d],
            rid = arr[1][1]? arr[1] : [arr[1]],
            sel = rid.map(function(x){ return '.gg' + arr[0] + ' .pp' + x});
        d3.selectAll(sel.join(',')).classed('magenta', false)
    })
}
function addCon(cons){
    cons.forEach(function(d){
        var arr = contigData[d],
            gid = arr[0],
            rid = arr[1][1]? arr[1] : [arr[1]];
        if (!$('.con'+d)[0]){
            var aa;
            if (rid.length==1){
                aa = acc_genome[gid][rid[0]]
            } else {
                rid.forEach(function(x){
                    var a = acc_genome[gid][x];
                    if (!aa){ aa = a}
                    else if (a.length>aa.length){ aa = a}
                })
            }
            drawSynteny(aa)
        }
        var sel = rid.map(function(x){ return '.gg' + gid + ' .pp' + x});
        d3.selectAll(sel.join(',')).classed('magenta', true)
    })
}
function drawSynteny(acc){
    var myStrain = d3.select("#allORF").append("g")
                .attr("class", 'con'+acc[0])
                .classed("orfs", true),
        line = myStrain.append("g").attr("class", "grayLine");

    var lTri=8, hRNA=2;
    var gid = contigData[acc[0]][0],
        n = gids.indexOf(gid),
        base = hUnit * n;        
    var start0 = 0;
    $.each(acc, function(i, ac){
        //draw base line
        var cc = contigData[ac];
        var pp1 = start0;
        if (cc[3]){ pp1 += cc[4]? cc[4] : gapCont }
        var pp2 = pp1 + cc[2],
            p1 = scaleSyn(pp1),
            p2 = scaleSyn(pp2);
        line.append("line").attr("x1", p1).attr("x2", p2).attr("y1", base).attr("y2", base);
        line.append("line").attr("x1", p1).attr("x2", p1).attr("y1", base-hORF*1.2).attr("y2", base+hORF*1.2);
        line.append("line").attr("x1", p2).attr("x2", p2).attr("y1", base-hORF*1.2).attr("y2", base+hORF*1.2);
            
        var oo = orfData.filter(function(d){ return d.con==ac});
            
        $.each(oo, function(z,obj){
            var end = Math.round(scaleSyn(obj.end*1 + pp1)*10)/10,
                start = Math.round(scaleSyn(obj.end*1-obj.L*1 + pp1)*10)/10;
            var xv = [start, end-lTri*(obj.L>0? 1 : -1), end, end-lTri*(obj.L>0? 1 : -1), start],
                hB = obj.rna? hRNA : hORF,
                yv = [base-hB/2, base-hB/2, base, base+hB/2, base+hB/2],
                points = $.map(xv, function(x,i){ return x + ',' + yv[i]}).join(' ');

            var tipText = obj.locus;
            if (obj.sym){ tipText += '&nbsp;&nbsp;[ <em>' + obj.sym + '</em> ]' }
            if (obj.rna){ tipText += '&nbsp;&nbsp;[ ' + obj.rna + 'RNA ]'}
               
            var co, lo;
            if (obj.cid && obj.cid>0 && obj.orth){ co = obj.cid + '-' + obj.orth }
            if (obj.cid && obj.cid<0){ lo = obj.locus.replace(/\./g, '')}
               
            var myclass =  obj.cid>0? (co && $.inArray(co, uvList)!=-1? 'ciduv' : cid_class[obj.cid]) : 'cidna';

            if (obj.cid){
                if (obj.cid<0){ myclass += ' ' + n+'_'+lo }
                else if (!obj.orth){ myclass += ' ' + obj.con + '_' + obj.cid}
            }

            var thisOrf = myStrain.append("polygon")
                    .attr("points", points)
				    .attr("class", myclass)
				    .on("mouseover", function() {
				        var xMouse = d3.event.pageX-45;
                        $("#replicontip").css("left", xMouse+"px").css("top", base+"px").show().html(tipText);
                        d3.select('#treeRep #s' + gid).classed("darkorange", true)
                    })
                    .on("mouseout", function() {
                        $("#replicontip").hide();
                        d3.selectAll('#treeRep text').classed("darkorange", false)
                    });
               
            if (!obj.orth) {
                thisOrf.classed('whiteFill',true);
                if (!obj.cid && !obj.rna) { thisOrf.classed('dotStroke',true)}
            }

            if (!obj.cid){ return }

            thisOrf.classed('finger',true)
                   .on("dblclick", function(){ getParalog(obj.cid>0? obj.cid:0, obj.locus) });
            
            if (obj.orth){
                thisOrf.attr("id", gid+'_'+co)
                        .attr('end', obj.end*1+pp1)
                        .attr('len', obj.L)
                        .on("click", function(){ accFromId(obj.cid, obj.orth, gid) })
            } else { thisOrf.on("click", function(){ locateORF(obj.cid>0? obj.cid : lo)}) }
        });
        start0 = pp2
    });
    myStrain.attr('len', start0)
}
function locateORF(cid){
    d3.selectAll('#synteny svg .shade:not(#shadeC)').remove();
    
    for (var i=0; i<currAcc.length; i++){
        var acc = currAcc[i],
            gid = contigData[acc][0],
            n = gids.indexOf(gid),
            base = hUnit*n,
            myorfs = $('#synteny .' + acc + '_'  + cid);
        
        if (!myorfs[0]){ continue }
        
        $.each(myorfs, function(i,myorf){
            var points = myorf.attributes.points.value.split(' ');
            var p1 = points[0].split(',')[0],
                p2 = points[2].split(',')[0],
                px = Math.min(p1,p2);
            
            px += orfTranslate[gid]? orfTranslate[gid] : 0;
            
            var shd = d3.select('#synteny svg>g').append("rect")
                .attr("class", 'shade dotStroke')
                .attr("x", px).attr("y", base-7)
                .attr("width", Math.abs(p2-p1)).attr("height", 14);
            
            if (!stdOid || stdOid && cid!=stdOid.split('-')[0]){
                shd.style("stroke", "gray")
            }
        })
    }
}

function queryLoc(query, exact){
    type =0;
    query = query.replace(/ /g, '');
    if (query.length == 0) {return}
    
    var url= myUrl + 'orf6/select?q=locus' + (exact? '' : 'F') + ':' + query + '&start=0&rows=' + (exact? 1 : 50) + '&wt=json&callback=?&json.wrf=' + (exact? 'getRs' : 'showSugg') + '&fl=cid,locus' + (exact? ',rna' : '');
    $.getJSON(url)
}
function querySym(query, exact){
    type=1;
    query = query.replace(/ /g, '');
    if (query.length == 0) {return}
    
    var url= myUrl + 'sym_cid6/select?q=sym' + (exact? '' : 'F') + ':' + query + '&start=0&rows=' + (exact? 1 : 10) + '&wt=json&callback=?&json.wrf=' + (exact? 'getRs' : 'showSugg');
    if (!exact){ url += '&fl=cid' + (exact? '' : ',sym')}
    $.getJSON(url)
}
function queryPfam(query){
    type=2;
    query = query.replace(/ /g, '');
    if (query.length == 0) {return}
    
    var url= myUrl + 'cdhit6/select?q=pfam:' + query + '&start=0&rows=13&wt=json&callback=?&json.wrf=getRs';
    $.getJSON(url)
}

function getRs(data){
    var docs = data.response.docs;
    var inf, obj;
    if (!docs.length){
        inf = 'not found in database'
    } else if (type==2 && docs.length>1){
        pfamId = $('#searchPfam').val();
        showSugg(docs); return
    } else {
        obj = docs[0];
        if (!type && !obj.cid){
            inf = 'is ' + (obj.rna? (obj.rna+'RNA') : 'pseudogene')
        }
    }
    
    if (inf){
        $('#newsPara')
            .html('<b>'+ (type==1? $('#searchSym').val() : (type==2? 'pf' + $('#searchPfam').val() : $('#searchLoc').val())) +'</b> '+ inf)
            .css("left", (type==1? 533 : (type==2? 810 : 210)) + 'px').show();
        return
    }
    getParalog(obj.cid>0? obj.cid : 0, type? '' : obj.locus)
}
function showSugg(data) {
    var cont = '';
    var docs = type==2? data : data.response.docs;
    $.each(docs, function(i, obj) {
        if (!type && !obj.cid){ return }
        cont += '<p onclick=\"getParalog(' + (obj.cid>0? obj.cid : 0) + (type? '' : ',\''+obj.locus+'\'') + ')\">' + (type==2? 'pf' + pfamId + '-' + (i+1) : (type==1? obj.sym : obj.locus)) + '</p>'
    });
    $('#suggestion').html(cont).show()
        .css("left", (type==2? 850 : (type? 533 : 210)) + "px")
        .css("height", type? "auto" : "400px")
        .css("overflow-y", type? "hidden" : "scroll")
}

var orfPara;
function getParalog(cid, myLocus) {
    $('#main_tabs').tabs({ active:2 });

    if (myLocus && locus && locus==myLocus){
        $('#searchSym').val(symbol);
        $('#searchLoc').val(locus);
        $('#searchPfam').val(pfamId);
        return
    }
    if (myLocus){locus = myLocus}
    
    if (cid && cdhit && cdhit==cid){
        if (myLocus){ changeHL() }
        $('#searchSym').val(symbol);
        $('#searchLoc').val(locus);
        $('#searchPfam').val(pfamId);
        return
    }
    
    $('#paraInf').show();
    var url = myUrl + 'orf6/select?q=';
    if (cid){
        cdhit = cid;
        url += 'cid:'+cdhit+'&start=0&rows=650'
    } else {
        cdhit = 0;
        url += 'locus:'+myLocus+'&start=0&rows=1' 
    }
    url += '&wt=json&json.wrf=?';

    $.getJSON(url, function(data){
        orfPara = data.response.docs;
        
        $('#paraName tr').remove();
        $('#treeP, #seqDown').empty();
        $(':radio').prop('checked',false);
        $('#selectNtAa').hide();

       if (cdhit){
           $('#dbl_click').show();
           var ur = myUrl + 'cdhit6/select?q=cid:'+cdhit+'&start=0&rows=1&wt=json&json.wrf=?';
            $.getJSON(ur, function(d){
                var cdhitData = d.response.docs[0];
                symbol = cdhitData.sym? cdhitData.sym.join(', ') : '';
                $('#searchSym').val(symbol);

                if (cdhitData.dnd){
                    paraDnd = cdhitData.dnd;
                    $('#lgScale span, #lgScale svg').show()
                    draw_treeP()
                } else { noDnd() }

                if (!myLocus){ locus = idsP[0] }
                $('#searchLoc').val(locus);
                locus_old = locus.toUpperCase().replace(/\./g, '');
                
                pfamId = cdhitData.pfam? cdhitData.pfam : '';
                $('#searchPfam').val(pfamId);

                paraPfam(cdhitData.pfam);
                if (cdhitData.patric) { paraPatric(cdhitData.patric) } else { $('#patric').hide() }
                if (cdhitData.pmed) { paraPub(cdhitData.pmed) } else { $('#pubmed').hide() }
                if (cdhitData.tigr) { paraTigr(cdhitData.tigr) } else { $('#tigr').hide() }
                paraName();
                setTabWid()
            })
        } else {
           $('#dbl_click').hide();
            noDnd();
            locus_old = locus.toUpperCase().replace(/\./g, '');
            paraName();
            setTabWid();
            $('#pfam, #patric, #tigr, #pubmed').hide()
        }

        $('#searchLoc').val(locus)
    });

    var w_limit = wTreeP - 4;
    function draw_treeP(){
        var root = d3.hierarchy(parseNewick(paraDnd), function(d) {return d.children})
            .sort(function(a,b){return d3.descending(a.data.length, b.data.length)});
        idsP = root.leaves().sort(function(a,b) { return (a.x - b.x)}).map(function(d){return d.data.name});
        orfPara = orfPara.sort(function(a,b){return idsP.indexOf(a.locus) - idsP.indexOf(b.locus)});
        hsvgP = idsP.length*spH;

        var cluster = d3.cluster().size([hsvgP, 0]).separation(function(a, b) { return 1 });
        cluster(root);
        resetX(root, root.data.length=0, w_limit/maxLength(root));

        var leftBound = 1;
        var svg = d3.select('#treeP').append("svg")
                .attr("width", leftBound + w_limit)
                .attr("height", hsvgP),
            chart = svg.append("g")
                .attr("transform", "translate(" + leftBound + "," + (-2) + ")");

        chart.append("g")
            .attr("class", "blackLine")
            .selectAll("path").data(root.links()).enter()
            .append("path")
            .attr("d", linkVariable);

        chart.append("g")
            .attr("class", "link-extensions")
            .selectAll("path")
            .data(root.links().filter(function(d) { return !d.target.children })).enter()
            .append("path")
            .attr("d", linkExtensionVariable);

    //draw scale:
        var unit = 30;
        if (!d3.select('#txtScale').html()){
            var lineScale = d3.select("#lgScale svg").attr("width",unit)
                    .append("g").attr("transform", "translate(3,6)");
            lineScale.append("line").attr("x1", 0).attr("x2", unit).attr("y1", 0).attr("y2", 0);
            lineScale.append("line").attr("x1", 0).attr("x2", 0).attr("y1", -2).attr("y2", 2);
            lineScale.append("line").attr("x1", unit).attr("x2", unit).attr("y1", -2).attr("y2", 2)
        }
        $("#txtScale").html(d3.format(".2g")(unit*maxLength(root)/w_limit))

        function linkExtensionVariable(d) {
            return linkStep(d.target.y, d.target.x, w_limit, d.target.x)
        }
    }

    function changeHL(){
        var lu = locus.toUpperCase().replace(/\./g, '');
        if (lu!==locus_old){
            $('#' + locus_old).css("background-color", "inherit");
            $('#' + lu).css("background-color", color_paraName);
            locus_old = lu
        }
    }
    function paraName(){
        var trs;
        $.each(orfPara, function(i,obj){
            var orth = obj.orth? obj.orth : 0,
                id = obj.locus,
                gid = contigData[obj.con][0],
                grpid = contigData[obj.con][1];
            if (!grpid[0]){ grpid = [grpid] }

            var lu = id.toUpperCase().replace(/\./g, '');
            trs += '<tr height="' + spH + 'px" id="' + lu + '"' + (orth? (' class="orth' + orth + '"') : '') + (lu==locus_old? (' style="background-color:' + color_paraName + '"') : '') + '><td>';
            if (orth){
                trs += '<a href="#" ondblclick="accFromId(' + cdhit  + ',' + orth + ')">'
            }
            trs += id;
            if (orth){ trs += '</a>' }
            trs += '</td><td>';

            var strainObj = strainData[gid];
            var sp = strainObj.spId;
            var ar = [strainObj.name, spcData[sp][0]];
            ar.push($.map(grpid, function(c){ return repName[c][0] }).join('|'));
		  trs += ar.join('</td><td>');
		  trs += '</td></tr>'
        });
        $('#paraName').append(trs);
        $('#paraName a').click(function(e){ e.preventDefault() })
    }
    function paraPfam(pfam){
        $('#pfamId').html(pfam? pfam : 'unclassified');
        var pfamInf;
        if (pfam){ pfamInf=famData.pfam[pfam] }
        $('#pfamInf').html(pfamInf? pfamInf : '')
        $('#pfam').show()
    }
    function paraPatric(patric){
        var trs;
        $.each(patric, function(i,id){
            var anno = famData.patric[id];
			trs += '<tr><td>' + id + '</td><td><div class="vr"></div></td><td>' + (anno? anno : '') + '</td></tr>'
		});
	
	   $('#patric tr:not(:first-child)').remove();
	   $('#patric').append(trs).show()
    }
    function paraPub(pub){
        var trs;
        $.each(pub, function(n,id){
            var obj = pubData[id];
            trs += '<tr><td>&bull;</td><td>' + (obj[0]? obj[0] + ' <em>et al</em>. ' : '[No authors listed]. ') + obj[3] + '. ' + obj[1] + ' <em>' + obj[2] + '</em>. <a href="' + dbLink.pubmed + id + '" target="_blank">PMID' + id + '</a></td></tr>'
        });
	   $('#pubmed tr:not(:first-child)').remove();
	   $('#pubmed').append(trs).show()
    }
    function paraTigr(tigr){
        var trs;
        var n=0;
        $.each(tigr, function(i,id){
            var obj = famData.tigr[id];
            if (n){ trs += '<tr><td colspan="2"><hr></td></tr>' }
			trs += '<tr><td>TIGR ID</td><td>' + id + '</td></tr>';
			trs += '<tr><td>Description</td><td>' + obj[0] + '</td></tr>';
			if (obj[1]){ trs += '<tr><td>Comment</td><td>' + obj[1] + '</td></tr>' }
			if (obj[2]){ trs += '<tr><td>DB Ref</td><td>' + obj[2].join(' | ') + '</td></tr>' }
			if (obj[3]){ trs += '<tr><td>Pub Ref</td><td>' + obj[3].join(' | ') + '</td></tr>' }
			if (obj[4]){ trs += '<tr><td>Role</td><td>' + $.map(obj[4], function(v,i){
                    var myRole = famData.role[v];
                    return (myRole[0]? ('Main: ' + famData.mainRole[myRole[0]] + '; ') : '') + 'Sub: ' + myRole[1]
                }).join('<br>') + '</td></tr>'
            }
			if (obj[5]){ trs += '<tr><td>GO</td><td>' + $.map(obj[5], function(v,i){return '<a href="' + dbLink.go + v + '" target="_blank">GO:'+v + '</a>'}).join(' | ') + '</td></tr>' }
            n++
		});
	
	   $('#tigr tr:not(:first-child)').remove();
	   $('#tigr').append(trs).show()
    }
    function noDnd(){
        idsP = orfPara.map(function(o){ return o.locus });
        hsvgP = idsP.length*spH;
        $('#lgScale span, #lgScale svg').hide()
    }
    function setTabWid(){
        $('#sub_tabs').tabs({ active:0 }).css('min-height', (hsvgP<457? 500 : hsvgP+43)+'px');
        var wSubTabP = 1012-wTreeP-$('#nameP').width();
        $('#sub_tabs').css('width', wSubTabP+'px');
        $('#seqP').css("width", (wSubTabP-wThumbnail-18) + 'px');
        $('#seqDown').css("width", (wSubTabP-25) + 'px')
    }
}


var cid_seq, seqNt, stdSeq, pLength;
function thumbnail(){
    if ($('#seqP').is(":visible")){
        $('#selectNtAa').show()
    } else {
        $('#selectNtAa').hide()
    }
    if (cid_seq && cid_seq==cdhit) {
        return}
    cid_seq = cdhit;

    $('#seqP, #selectNtAa').hide();
    $("#thumbnail, #seqPara").empty();
    
    seqNt = orfPara.map(function(d){ return d.alnNt});
    stdSeq = seqNt[0].split('');
    pLength = stdSeq.length;
    ratioP = pLength/wThumbnail;
    var lowerYbound = 6;
    var svg = d3.select("#thumbnail") .append("svg") .attr("width", wThumbnail+1) .attr("height", lowerYbound+hsvgP+2*spH);

    var line = svg.append("g").attr("class", 'tnStroke');
    var seg = svg.append("g").attr("class", 'orfPara');
    
    $.each(seqNt, function(i,seq){
        var py = (spH-hORF)/2+spH*i+1;
        line.append("rect")
            .attr("x", 0)
            .attr("y", py)
            .attr("width", wThumbnail)
            .attr("height", hORF);
        
        var ss = [], good, blank, aa=seq.split('');
        for (var n=0; n<aa.length; n++){
            if (aa[n]!='-' && !good){ ss.push(n); good=1; blank=0 }
            else if (aa[n]=='-' && !blank) {
                good=0; blank=1;
                if (n) { ss.push(n-1) }
            }
        }
        for (var j=0; j<ss.length; j+=2){
            var s1 = ss[j],
                s2 = j+1==ss.length? aa.length : ss[j+1];
            seg.append("rect")
                .attr("x", s1/ratioP)
                .attr("y", py)
                .attr("width", (s2-s1+1)/ratioP)
                .attr("height", hORF)
        }
    });
    
    var base = hsvgP+spH;
    line.append("line").attr("x1", 0).attr("x2", wThumbnail).attr("y1", base).attr("y2", base);
    line.append("line").attr("x1", 0).attr("x2", 0).attr("y1", base-2).attr("y2", base+2);
    line.append("line").attr("x1", wThumbnail).attr("x2", wThumbnail).attr("y1", base-2).attr("y2", base+2);

    svg.append("text")
        .attr("x",wThumbnail/2)
        .attr("y", base+12)
        .text(numberWithCommas(pLength)+' bp');

//draw guideline
    guidelineP = svg.append("line").attr("class","guideline")
    			.attr("y1",0).attr("y2",base)
    			.attr("x1","-50px").attr("x2","-50px")
}


var coden = {AAA:"K",AAG:"K",AAC:"N",AAT:"N",GAT:"D",GAC:"D",GAG:"E",GAA:"E",CAA:"Q",CAG:"Q",CAC:"H",CAT:"H",TAT:"Y",TAC:"Y",TAG:"*",TAA:"*",TGA:"*",TGG:"W",TGC:"C",TGT:"C",CGT:"R",CGC:"R",CGG:"R",CGA:"R",GGA:"G",GGG:"G",GGC:"G",GGT:"G",AGT:"S",AGC:"S",AGG:"R",AGA:"R",ACA:"T",ACG:"T",ACC:"T",ACT:"T",GCT:"A",GCC:"A",GCG:"A",GCA:"A",CCA:"P",CCG:"P",CCC:"P",CCT:"P",TCT:"S",TCC:"S",TCG:"S",TCA:"S",TTA:"L",TTG:"L",TTC:"F",TTT:"F",CTT:"L",CTC:"L",CTG:"L",CTA:"L",GTA:"V",GTG:"V",GTC:"V",GTT:"V",ATT:"I",ATC:"I",ATG:"M",ATA:"I"},
    colorAA = {G:"green", S:"green", T:"green", Y:"green", C:"green", Q:"green", N:"green", K:"blue", R:"blue", H:"blue", D:"red", E:"red"};
function drawAlnSeq(){
    var svg_w = pLength*unitW,
        hAxis = 20,
        svg_h = spH*idsP.length + hAxis;

    var svg = d3.select('#seqPara').append("svg").attr("width", svg_w).attr("height", svg_h);

    var codenBg = svg.append("g").attr("class", 'codenBg');
    for (var i=0; i<pLength/6; i++){
        codenBg.append("rect").attr("x",i*6*unitW).attr("y",0).attr("width",unitW*3).attr("height",svg_h)
    }

    svg.append("g").attr("id", "seqNt")
        .selectAll("text").data(seqNt).enter().append("text")
        .attr("transform", function(d,i) { return "translate(0," + (i*spH) + ")" })
        .attr("dy", "1em")
        .html(function(d,i) {
            if (i){ d = matchSeq(d)}
            return d.replace(/-/g, '&nbsp;')
        });

    var axis=svg.append("g").attr("transform", "translate(-3," + (spH*idsP.length+2) + ")"),
        tick = axis.append("g").attr("class", "tick"),
    	mark = axis.append("g").attr("class", "mark");
    for (var i=1; i<pLength/30; i++){
        var px = i*30*unitW;
        tick.append("line").attr("x1",px).attr("x2",px).attr("y2", 6);
        mark.append("text").attr("x",px+3).attr("y",hAxis-2).text(i*30)
    }
    $("#seqP").show().scrollLeft(0);
    $("#selectNtAa").html("Show AA").fadeIn(2000);
    showAA=1
}
function matchSeq(d){return arr = d.split('').map(function(a,i){ return a==stdSeq[i] && a!='-'? '.' : a }).join('')}

function downSeq(i){
	var content='';
	if (i==2){ if (paraDnd){content = paraDnd} }
	else {
        for (var n=0; n<orfPara.length; n++){
            var obj = orfPara[n],
                loc = obj.locus,
                seq = obj.alnNt,
                strainObj = strainData[contigData[obj.con][0]];
            content += '>' + spcData[strainObj.spId][0] + '|' + strainObj.name + '|' + loc + '<br>' + (i? transAA(seq) : seq) + '<br>'
        }
	}
	$('#seqDown').html(content)
}

function hideSelectBtn(){$("#selectNtAa").hide()}

var stdAA;
function aaSeqTbl(){
    var aaSeq = seqNt.map(function(d){return transAA(d)});
    stdAA=aaSeq[0].split('');
    
    d3.select("#seqPara svg").append("g").attr("id", "seqAA").attr("class", "opc0")
        .selectAll("text").data(aaSeq).enter().append("text")
        .attr("transform", function(d,i) { return "translate(" + unitW + "," + (i*spH) + ")" })
        .attr("dy", "1em")
        .html(function(d,i) {
            var arr =  matchSeqAA(d,i)
            return arr.join('&nbsp;&nbsp;').replace(/-/g, '&nbsp;')
        })
}
function transAA(s){
    var arr = s.split('');
    var aa = '';
    for (var i=0; i<arr.length; i+=3){
        var c = arr[i] + arr[i+1] + arr[i+2];
        aa += c=='---'? '-' : coden[c.toUpperCase()]
    }
    return aa
}
function matchSeqAA(d,j){
    var arr = d.split('').map(function(a,i){
        if (j && a==stdAA[i] && a!='-'){ a = '.' }
        var span = '';
        if (colorAA[a]){span = '<tspan class="' + colorAA[a] + '">'}
        return span + a + (span? '</tspan>' : '')
    })
    return arr
}


function drawTree(){
    var wTree = wsvgTree - 62,
        hsvgTree = hUnit*gids.length,
        lowerXbound=1,
        gap=6;
    
    var cluster = d3.cluster().size([hsvgTree, 0]).separation(function(a, b) { return 1 });
    var root = root2;
    cluster(root);
    resetX(root, root.data.length=0, wTree/maxLength(root));

    var svg = d3.select('#treeRep').append("svg")
            .attr("width", wsvgTree)
            .attr("height", hsvgTree)
            .append("g").attr("transform", "translate(" + lowerXbound + ",0)");

    svg.append("g")
        .attr("class", "whiteStroke")
        .selectAll("line").data(gids).enter()
        .append("line")
        .attr("transform", function(d,i) { return "translate(" + (-4) + "," + (i*hUnit) + ")" })
        .attr("x1", wTree+8).attr("x2", wsvgTree);

    svg.append("g")
        .attr("class", "blackLine")
        .selectAll("path").data(root.links()).enter()
        .append("path")
        .attr("d", linkVariable);

    svg.append("g")
        .attr("transform", "translate(" + (wTree+gap) + "," + (hUnit-5) + ")")
        .selectAll("text").data(gids).enter()
        .append("text")
        .attr("transform", function(d,i) { return "translate(0," + (i*hUnit) + ")" })
        .attr("id", function(d){return 's'+d})
        .text(function(d) { return strainData[d].name })
}
function resetX(d, y0, k) {
    d.y = (y0 += d.data.length) * k;
    if (d.children) d.children.forEach(function(d) { resetX(d, y0, k) })
}
var lineFunction = d3.line().x(function(d){return d[0]}).y(function(d){return d[1]}).curve(d3.curveLinear);
function linkStep(x1, y1, x2, y2) {
    var arr = [[x1,y1]];
    if (y2!=y1){ arr.push([x1,y2]) }
    arr.push([x2,y2]);
    return lineFunction(arr)
}
function linkVariable(d) {
    return linkStep(d.source.y, d.source.x, d.target.y, d.target.x)
}


function openLink(site,ids){
    if (!ids){ return }
    var win = window.open(dbLink[site] + ids.join(','), '_blank');
    if (win) { win.focus() } else { alert('Please allow popups for this website') }
}

function drawLg(){
    var svg_h = 19, svg_w=232, base=8, recH=4, recW=17;

    var lTxt = ['pseudogene or sequence error', 'RNA']

    var svg = d3.select('#legORF').append("svg").attr("width", svg_w).attr("height", svg_h);
    var rect = svg.append("g").attr("class", "tnStroke"),
        txt_leg = svg.append("g").attr("class", "font10 ita");

    var xPos=7;
    for (var i=0; i<2; i++){
		rect.append("rect")
            .attr("x", xPos)
            .attr("y", i? base+1 : base)
            .attr("width", recW)
            .attr("height", i? recH-2 : recH)
            .classed("dotStroke", i? false : true);
 		txt_leg.append("text").attr("x", xPos+22).attr("y", base+5).text(lTxt[i]);
        xPos += recW + 159
    }
}


function parseNewick(a){for(var e=[],r={},s=a.split(/\s*(;|\(|\)|,|:)\s*/),t=0;t<s.length;t++){var n=s[t];switch(n){case"(":var c={};r.children=[c],e.push(r),r=c;break;case",":var c={};e[e.length-1].children.push(c),r=c;break;case")":r=e.pop();break;case":":break;default:var h=s[t-1];")"==h||"("==h||","==h?r.name=n:":"==h&&(r.length=parseFloat(n))}}return r}

function maxLength(d) {
  return d.data.length + (d.children ? d3.max(d.children, maxLength) : 0)
}
function setRadius(d, y0, k) {
  d.radius = (y0 += d.data.length) * k;
  if (d.children) d.children.forEach(function(d) { setRadius(d, y0, k) })
}

function linkFun(d) { return linkStepC(d.source.x, d.source.radius, d.target.x, d.target.radius) }
function linkConstant(d) { return linkStepC(d.source.x, d.source.y, d.target.x, d.target.y) }

function linkExtensionFun(d) { return linkStepC(d.target.x, d.target.radius, d.target.x, innerR) }
function linkExtensionConstant(d) { return linkStepC(d.target.x, d.target.y, d.target.x, innerR) }

function linkStepC(startAngle, startRadius, endAngle, endRadius) {
  var c0 = Math.cos(startAngle = (startAngle - 90) / 180 * Math.PI),
      s0 = Math.sin(startAngle),
      c1 = Math.cos(endAngle = (endAngle - 90) / 180 * Math.PI),
      s1 = Math.sin(endAngle);
  return "M" + startRadius * c0 + "," + startRadius * s0
      + (endAngle === startAngle ? "" : "A" + startRadius + "," + startRadius + " 0 0 " + (endAngle > startAngle ? 1 : 0) + " " + startRadius * c1 + "," + startRadius * s1)
      + "L" + endRadius * c1 + "," + endRadius * s1
}

function testIta(d){
    return /\./.test(d) || d=='Dermacentor'
}

var scaleSpace=19;
function drawScale(svg, ruler, ir, max, x, y, id){
    var scaleL = ruler*ir/max;
    
    var scale = svg.append("g").attr("id",id),
        lineScale = scale.append("g")
            .attr("class", "tick")
            .attr("transform", "translate(" + (x-scaleL/2) + "," + (y-scaleSpace) + ")");

    lineScale.append("line").attr("x2",scaleL);
    lineScale.append("line").attr("y1",-2).attr("y2",2);
    lineScale.append("line").attr("x1",scaleL).attr("x2",scaleL).attr("y1",-2).attr("y2",2);

    scale.append("text").attr("class","midAnchor font10 mark")
        .attr("x", x).attr("y",y).text(ruler + ' AA Sub/Site')
}

function drawAxis(tick_mark, ratio, l, unit, sub, k) {
    var tick = tick_mark.append("g").attr("class", "tick"),
    	mark = tick_mark.append("g").attr("class", "mark"),
        hTick=4;

    tick.append("line").attr("x1", 1/ratio).attr("x2", l/ratio).attr("y1", base_tick).attr("y2", base_tick);
    tick.append("line").attr("x1", 1/ratio).attr("x2", 1/ratio).attr("y1", base_tick-hTick/2).attr("y2", base_tick+hTick).style("stroke-width", "1.6px");
    mark.append("text").attr("x",1/ratio).attr("y", base_tick-5).text(1).style("text-anchor", "start");

    for (i=1; i<=parseInt(l/unit); i++){
    	var x = i*unit/ratio;
		tick.append("line").attr("x1", x).attr("x2", x).attr("y1", i%5? base_tick : base_tick-hTick/2).attr("y2", base_tick+(i%5? hTick/2 : hTick));
		if (i%sub==0) { mark.append("text").attr("x", x).attr("y", base_tick-5).text(i*unit/k + (k==1000? 'k' : (k==1000000? 'M' : ''))) }
    }
}

function numberWithCommas(x) { return x.toString().replace(/\B(?=(\d{3})+(?!\d))/g, ",") }

function blast(){
    $('#newsBlast').hide();
    var seq = $('#sequence').val();
    
    seq = seq.replace(/ /g,'').replace(/-/g,'');
    if (!seq) { return }
    var id, ss='';
    var lines = seq.split('\n');
    $.each(lines, function(i,line){
        if (!line){ return }
        if (!id){
            if (!/^>/.test(line)){ ss += line }
            id = 'query'
        } else {
            ss += line
        }
    });
    if (!ss){ return }

    if (!/^[ACDEFGHIKLMNPQRSTVWY]*$/i.test(ss)){ $('#newsBlast').html('The sequence contains non-standard codes.').show() }
    
    $.ajax({
		url: "cgi-bin/runblast.pl",
        dataType: "json",
        data: 'ev=' + $('#eValue').val() + '&gid=' + $('#strainBlast').val() + '&seq=' + ss,
		error: function(){ alert('AJAX error.') },
		success: function(data) { showBlast(data, ss.length) }
    })
}
function showBlast(data, lenQ){
    if (data.rna){ $('#newsBlast').html('This is RNA sequence.').show(); return }
    else if (data.no){ $('#blastRes').html('<h3>No hits found</h3>'); return }
    $('#blastRes').empty();

    for (var i=0; i<data.length; i++){
        var obj = data[i];
        var gapQ = (obj.seqQ[2].match(/-/g) || []).length,
            gapS = (obj.seqS[2].match(/-/g) || []).length;
        
        var hitLocus = '';
        if (!obj.orfs){ hitLocus = obj.locus }
        else {
            var id = obj.orfs.shift();
            var cont = repName[id][0];
            if (id==31){ cont = repName[25][0] + '|' + cont }
            
            $.each(obj.orfs, function(n,o){ hitLocus += (n? ' - ' : '') + (o[0]? o[0] : (n? 'end' : 'start')) })
        }
        
        var res = '<div>';
        res += '<p><span>Subject: <b>';
        if (!i && (!obj.orfs || obj.orfs.length==1)){ res += '<a href="#" ondblclick="queryLoc(\'' + hitLocus  + '\',1)" onMouseOver="showInfo()" onMouseOut="hideInfo()">'}
        res += hitLocus;
        if (!i){ res += '</a>'}
        res += '</b>';
        if (cont){ res += ' at <b>' + cont + '</b>' }
        res += '</span>';

        res += '<span>E value: ' + obj.ev + '</span>';
        res += '<span>Positives: ' + Math.round(obj.Positives*100) + '%</span>';
        res += '<span>Identities: ' + Math.round(obj.Identities*100) + '%</span>';
        res += '<span>Gaps in Query: ' + gapQ + '</span>';
        res += '<span>Gaps in Subject: ' + gapS + '</span></p>';

        res += '<table><tr>';
        res += '<td class="tickSeq"><br>Query<br><br>Subject</td>';
        
        res += '<td><div class="blastSeqArea">';
        res += '<span class="tickSeq">' + obj.seqQ[0] + Array(obj.seqQ[1]-obj.seqQ[0]+gapQ-(obj.seqQ[0].toString().length-1)).join("&nbsp;") + obj.seqQ[1] + '</span><br>';
        res += '<span class="blastSeq">' + obj.seqQ[2] + '<br>' + modiSeqH(obj.seqH) + '<br>' + obj.seqS[2] + '</span><br>';
        res += '<span class="tickSeq">' + obj.seqS[obj.rev?1:0] + Array(obj.seqS[1]-obj.seqS[0]+gapS-(obj.seqS[0].toString().length-1)).join("&nbsp;") + obj.seqS[obj.rev?0:1] + '</span><br><br>';
        res += '</div></td>';
        
        res += '<td class="blastThumb" id="bt_' + i + '"></td>';
        res += '</tr></table>';
        res += '</div>'
        $('#blastRes').append(res);

        drawBlastThumb(i, obj.lenS, obj.seqS[0], obj.seqS[1], lenQ, obj.seqQ[0], obj.seqQ[1], obj.rev, obj.locus, obj.orfs)
    }
}
function modiSeqH(seq){ return seq.replace(/[A-Y]/ig,'|').replace(/\s/g, '&nbsp;') }
function showInfo(){ $('#infoBlast').show() }
function hideInfo(){ $('#infoBlast').hide() }
function drawBlastThumb(i, lenS, s1, s2, lenQ, q1, q2, rev, hitLocus, orfs){
    var wBlast = 260;
    if (orfs){
        var s0 = orfs[0]? orfs[0][1] : 0;
        s1 -= s0;
        s2 -= s0;
        lenS = (orfs[orfs.length-1]? orfs[orfs.length-1][2] : lenS) - s0
    }
    var p0=Math.max(q1,s1),
        len = p0 + Math.max(lenQ-q1, lenS-s1),
        mar = 16;
        ratioB = len/(wBlast-mar*2),
        baseQ=23, baseS=56, hRect=5, lArr=16,
        yv = [baseS-hRect, baseS-hRect, baseS-1, baseS+1, baseS+hRect, baseS+hRect];
    
    var xscale = d3.scaleLinear().domain([0, len]).range([mar, wBlast-mar]),
        off = xscale(p0);

    var svg = d3.select("#bt_"+i).append("svg")
            .attr("width",wBlast).attr("height", 93);

    var unmatch = svg.append("g").attr("class", "blastLine blastUnmatch"),
        orfRect = svg.append("g").attr("class", "blastOrf"),
        match = svg.append("g").attr("class", "blastLine blastMatch"),
        guide = svg.append("g").attr("class", "blastUnmatch");
    
    if (!orfs){
        var start = off - s1/ratioB,
            end   = off + (lenS-s1)/ratioB;
        var xv = [start, end-lArr, end, end, end-lArr, start],
            points = $.map(xv, function(x,i){ return x + ',' + yv[i]}).join(' ');
        
        orfRect.append("polygon").attr("points", points);
        svg.append("text").attr("x", off+(lenS-s1*2)/ratioB/2).attr("y",72).text(hitLocus)
    } else {
        $.each(orfs, function(n,o){
            if (!o){ return }
            var start = off+(o[1]-s0-s1)/ratioB,
                end   = off+(o[2]-s0-s1)/ratioB;
            var xv = o[3]? [start, end-lArr, end, end, end-lArr, start] : [end, start+lArr, start, start, start+lArr, end],
                points = $.map(xv, function(x,i){ return x + ',' + yv[i]}).join(' ');

            orfRect.append("polygon").attr("points", points)
            svg.append("text").attr("x", start+(o[2]-o[1])/ratioB/2).attr("y",72).text(o[0])
        })
    }

    unmatch.append("line").attr("x1", off-q1/ratioB).attr("x2", off+(lenQ-q1)/ratioB).attr("y1", baseQ).attr("y2", baseQ);
    unmatch.append("line").attr("x1", off-s1/ratioB).attr("x2", off+(lenS-s1)/ratioB).attr("y1", baseS).attr("y2", baseS);
    
    guide.append("line").attr("x1", off).attr("x2", rev? (off+(s2-s1)/ratioB) : off).attr("y1", baseQ).attr("y2", baseS);
    guide.append("line").attr("x1", off+(q2-q1)/ratioB).attr("x2", rev? off : off+(s2-s1)/ratioB).attr("y1", baseQ).attr("y2", baseS);

    match.append("line").attr("x1", off).attr("x2", off+(q2-q1)/ratioB).attr("y1", baseQ).attr("y2", baseQ);
    
//draw axis
    var xAxis = d3.axisBottom().scale(xscale).ticks(5).tickSize(-3, 0)
                .tickFormat(function(f){ return f + (s0? s0 : 0) });
    svg.append("g").attr("class", "xaxis").attr("transform", "translate(" + (p0-s1)/ratioB + ",82)").call(xAxis)
}


var inR=285, osSvg, osCluster;
function ospcTree(){
    d3.tsv("js-css/ospc.tsv", function(data) { ospc=data });

    var lTxt = 50,
        outerR=inR + lTxt,
        wsvg = outerR*2;

    $.ajax({
        async:false, dataType:"text", url:"js-css/ospc.dnd",
        success: function(data){
            osRoot = d3.hierarchy(parseNewick(data), function(d){ return d.children })
                .sum(function(d){ return d.children ? 0 : 1 })
                .sort(function(a,b){return (a.value-b.value) || d3.ascending(a.data.length, b.data.length)})
        }
    });

    osRoot.data.length = 0;

    var i = 0;
    osRoot.descendants().forEach(function(d){
        d.id=++i;
        var name = d.data.name;
        if (name && /cut/.test(name)){ d.gp = name.replace(/cut_/, 'G') }
    });

//collapse when node has Id
    var node_coll = osRoot.descendants().filter(function(d){ return d.gp });
    node_coll.forEach(function(dd){
        dd._children = dd.children;
        dd.children = null
    });

    var svg = d3.select("#ospcTree").append("svg").attr("width", wsvg).attr("height", wsvg);
    osSvg = svg.append("g").attr("transform","translate("+outerR+","+outerR+")");
    osSvg.append("rect").attr("id","rectHL").attr("class","hidden");

    osCluster = d3.cluster().size([360, inR]).separation(function(a,b){return 1});
    updateOs(osRoot,1);
    $('#collapse').attr("disabled", "disabled");

//legend
    $('#lgSpc svg').clone().appendTo('#lgOspSpc');
    $('#lgArea svg').clone().appendTo('#lgOspArea');

    d3.select("#lgOspSpc svg").append("rect")
        .attr("id", "rect_sp")
        .attr("class", "hidden")
        .attr("height", hLg-4)
        .style("stroke", "darkgray")

    drawScale(svg, 0.05, inR, maxLength(osRoot), 43, wsvg-35, "ospScale")
}

function updateOs(source, start) {
    osCluster(osRoot);
    setRadius(osRoot, 0, inR/maxLength(osRoot));

    if (start){ osRoot.x0 = osRoot.x; osRoot.radius0 = osRoot.radius }

    var nodes = osRoot.descendants(),
        links = osRoot.links();

    var node = osSvg.selectAll('g.node').data(nodes, function(d){return d.id });

    var nodeEnter = node.enter().append('g')
        .attr('class', 'node finger')
        .attr("id", function(d){ return d.data.name && !d.gp? d.data.name : '' })
        .attr("transform", function(d){return "rotate(" + (source.x0-90) + ")translate(" + source.radius0 + ",0)"})
        .on('click', osClick)
        .on('mouseover', function(d){ospSpcTip(d.data.name)})
        .on('mouseout', function(){ospSpcTip(false)});

    nodeEnter.append('circle')
        .style("fill", function(d) { return d._children? "black" : "#fff" });

    nodeEnter.append('text')
        .attr("dy", ".35em")
        .attr("x", function(d){ return ((d.children || d._children) &&  d.x<180) || ((!d.children && !d._children) &&  d.x>=180) ? -6 : 6 })
        .attr("text-anchor", function(d) { return ((d.children || d._children) &&  d.x<180) || ((!d.children && !d._children) &&  d.x>=180)? "end" : "start" })
        .attr("transform", function(d) { return "rotate(" + (d.x<180? 0 : 180) + ")" })
        .text(function(d) {return d.data.name? (d.gp? d.gp : d.data.name) : '' })

    var nodeUpdate = nodeEnter.merge(node);

    var duration = 1000;
    nodeUpdate
      .transition().duration(duration)
      .attr("transform", function(d){return "rotate(" + (d.x-90) + ")translate(" + d.radius + ",0)"});

    nodeUpdate.select('circle')
        .attr('r', function(d){ return d.data.name || d==osRoot? 3 : 0})
        .style("fill", function(d) {
            var sp,
                condi = d.data.name && !d.gp;
            if (condi){
                var spp = ospSeq[d.data.name][0];
                if (spp){ sp = spcData[spp][1] }
            }
            return condi? (sp? color_sp(sp) : 'darkgray') : (d._children? "black" : "#fff") })
        .style("stroke", function(d) { return d.data.name && !d.gp? "none" : (d==osRoot? "darkgray" : "black") })
        .attr('cursor', function(d){ return d==osRoot? 'auto' : 'pointer' });

    var nodeExit = node.exit()
        .transition().duration(duration)
        .attr("transform", function(d){return "rotate(" + (source.x-90) + ")translate(" + source.radius + ",0)"})
        .remove();

    nodeExit.select('text').style('fill-opacity', 1e-6);

    var link = osSvg.selectAll('path.ospcLink')
            .data(links, function(d){ return d.target.id });

    var linkEnter = link.enter().insert('path', "g")
        .attr("class", "ospcLink")
        .attr('d', function(d){ return ospcStepX(source.x0, source.radius0, d.source.x<=d.target.x) })

    linkEnter.merge(link)
        .transition().duration(duration)
        .attr('d', function(d){
            if (start){ return ospcStep(d.source.x, d.source.radius, d.target.x, d.target.radius) }
            var wasL = $(this).attr('d').split(' ')[3]*1;
            var isL = d.target.x>=d.source.x? 1 : 0;
            if (wasL==isL){ return ospcStep(d.source.x, d.source.radius, d.target.x, d.target.radius) }
            return ospcStep(d.source.x, d.source.radius, d.target.x, d.target.radius, wasL? 'l' : 's')
        });

    var linkExit = link.exit()
            .transition().duration(duration)
            .attr('d', function(d){ return ospcStepX(source.x, source.radius, d.source.x0<=d.target.x0) })
            .remove();

    nodes.forEach(function(d){
        d.x0 = d.x;
        d.radius0 = d.radius
    });

    function ospcStepX(angle, radius, isL) {
        var c = Math.cos(angle = (angle-90)/180 * Math.PI),
            s = Math.sin(angle);
        return "M" + radius * c + "," + radius * s
                + ("A" + radius + "," + radius + " 0 0 " + (isL? 1 : 0) + " " + radius * c + "," + radius * s)
                + "L" + radius * c + "," + radius * s
    }

    function ospcStep(startAngle, startRadius, endAngle, endRadius, ls) {
        var c0 = Math.cos(startAngle = (startAngle-90)/180 * Math.PI),
            s0 = Math.sin(startAngle),
            c1 = Math.cos(endAngle = (endAngle-90)/180 * Math.PI),
            s1 = Math.sin(endAngle);

        var isL;
        if (ls){ isL = ls=='l'? 1 : 0 }
        else { isL = endAngle>=startAngle ? 1 : 0 }

        return "M" + startRadius * c0 + "," + startRadius * s0
            + "A" + startRadius + "," + startRadius + " 0 0 "+isL+" " + startRadius*c1 + "," + startRadius*s1
            + "L" + endRadius * c1 + "," + endRadius * s1
        }
}
function osClick(d){
    $('#searchAllele, #searchIsolate').val('');
    if (d==osRoot){ return }
    if (d.data.name && !d.gp) { osInf(d.data.name); return }
    
    $('#ospInf').slideUp("slow")
    d3.selectAll('#rect_sp, #rectHL').classed("hidden", true);
    $("#ospcTree text:not(.mark)").css("font-weight", "normal").css("font-size","11px");
    toggle(d)
}
function toggle(d){
    if (d.children) {
        d._children = d.children;
        d.children = null
    } else {
        expand(d);
        d._children = null;
        $('#collapse').removeAttr("disabled")
    }
    updateOs(d)
}
function expand(d) {
    if (d._children){
        d.children = d._children;
        d.children.forEach(expand)
    }
}
function collapse(){
    osRoot.descendants().filter(function(d){return d.data.name && d.children}).forEach(toggle);
    d3.selectAll('#ospcTree g.node, #ospcTree path.ospcLink').remove();
    updateOs(osRoot,1);
    $('#ospInf').slideUp("slow");
    d3.selectAll('#rect_sp, #rectHL').classed("hidden", true)
    $('#searchIsolate, #searchAllele').val('');
    $('#collapse').attr("disabled", "disabled")
}
function searchI(){
    $('#searchAllele').val('');
    var al = $('#searchIsolate').val();
    var os = ospc.filter(function(o){ return o.strain && o.strain.toUpperCase()==al.toUpperCase()});
    if (!os.length){ alert('No this isolate'); return }
    aaid2show(os[0].seqId)
}

function searchA(){
    $('#searchIsolate').val('');
    var al = $('#searchAllele').val();
    var os = ospc.filter(function(o){ return o.allele && o.allele.toUpperCase()==al.toUpperCase()});
    if (!os.length){ alert('No this allel name'); return }
    aaid2show(os[0].seqId)
}

function aaid2show(aaid){
    var dd = osRoot.descendants().filter(function(d){ return d.data.name==aaid });

    if (dd[0]){ osInf(aaid); return }

    var gp = ospSeq[aaid][1];
    dd = osRoot.descendants().filter(function(d){ return d.gp=='G'+gp });
    toggle(dd[0]);
    osInf(aaid,1)
}

function osInf(aaId, norect){
    $("#ospcTree text:not(.mark)").css("font-weight", "normal").css("font-size","11px").css("fill", "#b3b3b3");
    $('#'+aaId+' text').css("font-weight", "bolder").css("font-size","14px").css("fill", "black");

    if (norect){
        d3.select('#rectHL').classed("hidden",true)
    } else {
        var sel = d3.select('[id="' + aaId + '"]');
        var bbox = sel.select('text').node().getBBox();
        d3.select('#rectHL')
            .attr("width", bbox.width+4)
            .attr("height", bbox.height+2)
            .attr("transform", sel.attr("transform") + "translate(4,-9)")
            .classed("hidden",false)        
    }

    $('#ospName').html(aaId);

    var sp = ospSeq[aaId][0];
    if (sp){
        var pos = $.inArray(spcData[sp][1], spcGrpList);
        d3.select("#rect_sp").attr("y", hLg*pos+1).classed("hidden", false)
    } else {
        d3.select("#rect_sp").classed("hidden", true)
    }
    
    var os = ospc.filter(function(d){ return d.seqId==aaId});
    var trs = '';
    os.forEach(function(o){
        var istmp = /tmp/.test(o.acc);
        trs += '<tr><td>' + (istmp? '' : '<a href="' + dbLink.ncbi + o.acc + '" target="_blank">') + "&#9643" + (istmp? '' : '</a>') + '</td>';
        trs += '<td class="geo' + (o.geo? o.geo : 0) + '" style="color:' + (o.geo? color_area(geoData[o.geo].areaId) : 'darkgray') + '">' + "&#9642" + '</td>';
        var source = o.vector? bioSource[o.vector] : (o.host? bioSource[o.host] : ''),
            isIta = testIta(source);
        trs += '<td' + (isIta? ' class="ita" ' : '') + '>' + source + '</td>';
        trs += '<td>' + o.strain + '</td><td>' + o.allele + '</td>';
        trs += '<td>' + (o.pmed? '<a href="' + dbLink.pubmed + o.pmed + '" target="_blank">' : '') + (o.pmed? "&#9643" : '') + (o.pmed? '</a>' : '') + '</td></tr>'
    });
    $('#ospInf table tr:not(:first-child)').remove();
    $('#ospInf table').append(trs);
    $("#ospInf").slideDown("slow");
    $('#ospInf td:nth-of-type(2)')
        .on('mouseover', function(e){ospGeoTip(e, this)})
        .on('mouseout', function(){ospGeoTip(false)})
}

function ospGeoTip(e, active){
    if (!active){$("#ospGeotip").hide(); return}
    var geoId = active.className.substring(3) * 1;
    var xMouse = e.pageX-9,
        yMouse = e.pageY - 100;
    $("#ospGeotip")
        .html(geoId? geoData[geoId].name : 'unknown')
        .css("left", xMouse+"px").css("top", yMouse + "px").show()    
    return;
    var obj = strainData[active],
        sp = obj.spId,
        isIta = testIta(obj.source);

    $("#tipSp").html(spc[sp][0]);
    $("#tipSo").html(obj.source? ((isIta? '<i>' : '') + obj.source + (isIta? '</i>' : '')) : '');
    $("#tipGe").html(geoData[obj.geoId].name);
    if (obj.cit){
        $("#tipCt").html(obj.cit.map(function(pid){
            var o = pubData[pid];
            return o[0]+ ' et al. (' + o[3] + ') <i>' + o[2] + '</i>'
        }).join("; &nbsp;"))
    } else { $("#tipCt").empty() }
    $("#straintip").css("left", xMouse+"px").css("top", yMouse + "px").show()    
}
function ospSpcTip(aaId){
    if (!aaId || /cut/.test(aaId)){$("#ospSpctip").hide(); return}
    var spid = ospSeq[aaId][0],
        xMouse = d3.event.pageX-13,
        yMouse = d3.event.pageY-110;
    $("#ospSpctip")
        .html('species: ' + (spcData[spid]? ('<i>'+spcData[spid][0]+'</i>') : 'Unnamed'))
        .css("left", xMouse+"px").css("top", yMouse + "px").show()    
}
function ospBlast(){
    $('#newsOspBlast').hide();
    var seq = $('#ospSeq').val();
    seq = seq.replace(/ /g,'').replace(/-/g,'');
    if (!seq) { return }
    var id, ss='';
    var lines = seq.split('\n');
    $.each(lines, function(i,line){
        if (!line){ return }
        if (!id){
            if (!/^>/.test(line)){ ss += line }
            id = 'query'
        } else {
            ss += line
        }
    });
    if (!ss){ return }

    if (!/^[ACDEFGHIKLMNPQRSTVWY]*$/i.test(ss)){ ospBlastNews('The sequence contains non-standard codes.'); return }
    
    $.ajax({
		url: "cgi-bin/ospBlast.pl",
        dataType: "json",
        data: 'seq=' + ss,
		error: function(){ alert('AJAX error.') },
		success: function(data) {
            if (data.error){ alert(data.error); return }
            if (data.rna){ ospBlastNews('This is RNA sequence.'); return }
            if (data.no){ ospBlastNews('No hits found'); return }
            var aaid = data.id;
            if (data.Gaps || data.Identities!=1 || data.Positives!=1){
                ospBlastNews('The closest match is ' + aaid + '.'); return
            }
            aaid2show(aaid)
        }
    })
}
function ospBlastNews(inf){
    $('#newsOspBlast').html(inf).show();
    $('#ospInf').slideUp("slow")
}


function loadRepliconMenu(){
    $.each(groupT, function(id,obj){
        var name = repName[id][0];
        if (id==32){ name = repName[26][0] + '|' + name }
        $('#repliconMenu').append('<option value="'+id+'"' + (id==2? '  selected="selected"' : '') + '>'+ name + '</option>')
    })
}
function addTransTitle(){
    var transTitle0 = [['c-di-GMP',25987708], ['Growth Phases',27706236,3], ['Fitness in Carbohydrates',27279039,4], ['Life stages',25425211,3]];
    var tdAdd;
    for (var i=0; i<transTitle0.length; i++){
        tdAdd += i? '<td class="spacer"></td>' : 0;
        tdAdd += '<th' + (transTitle0[i][2]? ' colspan="' + transTitle0[i][2] + '"' : '') + '><a href="' + dbLink.pubmed + transTitle0[i][1] + '" target="_blank">' + transTitle0[i][0] + '</a></th>'
    }
    $('#transTitle0').append(tdAdd)
}
function addTransContainer(){
    var trAdd = '<tr><td id="transLocus"></td>';
    for (var i=0; i<transTitle.length; i++){
        trAdd += '<td id="' + transTitle[i] + '" class="transBar"></td>';
        if (!i || i==3 || i==7){ trAdd += '<td></td>' }
    }
    trAdd += '</tr>';
    $('#svgTrans').append(trAdd);
    
    var div = $('<div>').attr('id', 'guidelineT_area');
    $('#svgTrans').append(div)
}
function draw_transbar(c) {
    var transOrf = groupT[c];

    $('#transLocus').empty();
    $.each(transOrf, function(n,l){
        $('#transLocus').append('<p><a href="#" ondblclick="queryLoc(\'' + l  + '\',1)">' + l + '</a></p>')        
    });
        
    var svg_h = topEdge + hUnitT * (transOrf.length-1) + 6;
    var y1 = topEdge-6, y2 = svg_h;

    for (var i=0; i<transTitle.length; i++){
        var infl = transTitle[i],
            inflRange = infl.replace(/[0-9]/g, ''),
            rg = rangeT[inflRange];
        
        $('#'+transTitle[i]).empty();

        var svg = d3.select('#' + transTitle[i]).append("svg").attr("width", wBar).attr("height", svg_h);

        var scale = d3.scaleLinear().domain([-rg,rg]).range([1, wBar-2]);
        
        svg.append("rect").attr("class", "bg_trans")
	       .attr("x", scale(-rg))
           .attr("y", y1)
           .attr("width",scale(rg*2))
           .attr("height", y2-y1);
    
        //draw axis
        var xAxis = d3.axisTop().scale(scale).ticks(5).tickSize(-svg_h, 0);
        svg.append("g").attr("class", "xaxis")
            .attr("transform", "translate(0,18)")
            .call(xAxis);

        svg.append("line").attr("class", 'zeroline')
            .attr("x1", scale(0)).attr("x2", scale(0))
            .attr("y1", y1).attr("y2", y2);

        $.each(transOrf, function(n, l) {
            if (!transData[l] || !transData[l][infl]){ return }
            var vv = transData[l][infl];
            var base = topEdge + hUnitT * n;
        
            var v = i<8? vv[0] : vv;
            svg.append("line")
                .attr("x1", scale(0)).attr("x2", scale(v))
                .attr("y1", base).attr("y2", base)
                .attr("class", vv[1] || i>7? "transbarSig" : "transbarNo")
        })
    }
    
    //draw guideline
    $("#guidelineT_area").empty();
    var w = $('#svgTrans').width()? $('#svgTrans').width() : 1020;
    var svg = d3.select("#guidelineT_area")
              .append("svg")
              .attr("width", w-85)
              .attr("height", ($('#svgTrans').height()? $('#svgTrans').height() : 532)-62);
    guidelineT = svg.append("line").attr("class","guideline")
    			.attr("x1",0).attr("x2",w)
    			.attr("y1","-50px").attr("y2","-50px")
}
function drawSigLeg(){
    var svg_h = 19, svg_w=196, base=9, recH=4, recW=20;

    var recCls=['transbarSig', 'transbarNo'],
        lTxt = ['signifiicant', 'non-significant']

    var svg = d3.select('#legSig').append("svg").attr("width", svg_w).attr("height", svg_h);
    var txt_leg = svg.append("g").attr("class", "font10 ita");
    
    var xPos=7;
    for (var i=0; i<2; i++){
		svg.append("line")
            .attr("x1", xPos).attr("x2", xPos+recW)
            .attr("y1", base).attr("y2", base)
            .attr("class", recCls[i]);
 		txt_leg.append("text").attr("x", xPos+24).attr("y", base+3).text(lTxt[i]);
        xPos += recW + 73
    }
}

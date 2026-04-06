function _1(md){return(
md`# Buble world map`
)}

function _map_orign(htl,width,height,path,outline,graticule,land,borders){return(
htl.html`<svg viewBox="0 0 ${width} ${height}" style="display: block;">
    <path d="${path(outline)}" fill="#fff"></path>
    <path d="${path(graticule)}" stroke="#eaf4fc" fill="none"></path>
    <path d="${path(land)}"></path>
    <path d="${path(borders)}" fill="none" stroke="#FFF"></path>
    <path d="${path(outline)}" fill="none" stroke="#000"></path>
</svg>`
)}

function _borders(topojson,world){return(
topojson.mesh(world, world.objects.countries, (a, b) => a !== b)
)}

function _projection(d3){return(
d3.geoEqualEarth()
)}

function _path(d3,projection){return(
d3.geoPath()
    .projection(projection)
)}

function _outline(){return(
{type: "Sphere"}
)}

function _height(d3,projection,width,outline)
{
  const [[x0, y0], [x1, y1]] = d3.geoPath(projection.fitWidth(width, outline)).bounds(outline);
  const dy = Math.ceil(y1 - y0), l = Math.min(Math.ceil(x1 - x0), dy);
  projection.scale(projection.scale() * (l - 1) / l).precision(0.2);
  return dy;
}


function _graticule(d3){return(
d3.geoGraticule10()
)}

function _land(topojson,world){return(
topojson.feature(world, world.objects.land)
)}

function _svg(d3,map_orign){return(
d3.select(map_orign)
)}

function _world(FileAttachment){return(
FileAttachment("countries-50m.json").json()
)}

function _12(svg,topojson,world,path){return(
svg.append('path')
    .datum(topojson.feature(world, world.objects.land))
    .attr('fill', '#828384')
    .attr('d', path)
)}

function _13(svg,topojson,world,path){return(
svg.append('path')
    .datum(topojson.mesh(world, world.objects.countries))
    .attr('fill', 'none')
    .attr('stroke', 'none')
    .attr('stroke-linejoin', 'round')
    .attr('d', path)
)}

function _map(__query,FileAttachment,invalidation){return(
__query(FileAttachment("map_addice@1.csv"),{from:{table:"map_addice"},sort:[],names:[{name:"Pfam",column:"pfam"}],slice:{to:null,from:null},filter:[],select:{columns:["ID","pfam","Longitude","Latitude","Habitat","Hit"]}},invalidation)
)}

function _radius(d3,map){return(
d3.scaleSqrt()
    .domain([0, d3.max(map, d => d.Hit)])
    .range([5, 15])
)}

function _color(d3,map){return(
d3.scaleOrdinal()
    .domain([...new Set(map.map(d => d.Habitat))])
    .range(d3.schemeCategory10)
)}

function _color_1(d3){return(
d3.scaleOrdinal()
.domain([
  "Bioreactor",
  "Solid",
  "Wastewater",
  "Other engineered",
  "Freshwater",
  "Marine",
  "Thermal",
  "Soil",
  "Other enviroment",
  "Digestive",
  "Plants",
  "Microbial",
  "Other Host-associated",
  "Unclassfied"
])
.range([
  "#8dd3c7",
  "#ffffb3",
  "#bebada",
  "#D1D1B1",
  "#b3de69",
  "#80b1d3",
  "#fdb462",
  "#FFB6C1",
  "#ece8d9",
  "#fb8072",
  "#ccebc5",
  "#bc80bd",
  "#d4d6c8",
  "#808080"
])
)}

function _bubbles(svg,map,projection){return(
svg.append('g')
    .selectAll('g')
    .data(map)
    .join('g')
      .attr('transform', d => `translate(${projection([d.Longitude, d.Latitude])})`)
)}

function _19(bubbles,d3,radius,color_1){return(
bubbles.each(function(d) {
  const bubble = d3.select(this)
  const r = radius(d.Hit)
  if (d.Pfam === 'PF00905') {
    bubble.append('circle')
      .attr('r', r)
      .attr('fill', color_1(d.Habitat))
      .attr('fill-opacity', 0.5)
  } else if (d.Pfam === 'PF00144') {
    bubble.append('polygon')
      .attr('points', `${-r}, ${r} ${r}, ${r} 0, ${-r}`)
      .attr('fill', color_1(d.Habitat))
      .attr('fill-opacity', 0.5)
  } else if (d.Pfam === 'PF13354') {
    bubble.append('rect')
      .attr('x', -r)
      .attr('y', -r)
      .attr('width', r * 2)
      .attr('height', r * 2)
      .attr('fill', color_1(d.Habitat))
      .attr('fill-opacity', 0.5)
  }
})
)}

function _legendHabitat(svg){return(
svg.append("g")
  .attr("transform", "translate(0, 0)")
)}

function _21(legendHabitat,color_1){return(
legendHabitat.selectAll("rect")
  	.data(color_1.domain())
  	.enter().append("rect")
    	.attr("x", 10)
    	.attr("y", (d, i) => i * 20)
    	.attr("width", 20)
    	.attr("height", 10)
    	.attr("fill", color_1)
)}

function _22(legendHabitat,color_1){return(
legendHabitat.selectAll("text")
  	.data(color_1.domain())
  	.enter().append("text")
    	.attr("x", 30)
    	.attr("y", (d, i) => i * 20 + 6.5)
    	.attr("dy", "0.35em")
    	.text(d => d)
)}

function _legendPfam(svg){return(
svg.append('g')
  .attr("transform", "translate(30, 300)")
)}

function _pfamShapes(){return(
[
  {pfam: 'PF00905', shape: 'circle'},
  {pfam: 'PF00144', shape: 'polygon'},
  {pfam: 'PF13354', shape: 'rect'}]
)}

function _25(legendPfam,pfamShapes,d3){return(
legendPfam.selectAll('g')
    .data(pfamShapes)
    .join('g')
      .attr('transform', (d, i) => `translate(0, ${i * 20})`)
      .each(function(d) {
        const legend = d3.select(this)
        if (d.shape === 'circle') {
          legend.append('circle')
            .attr('cx', -9.5)
            .attr('cy', 9.5)
            .attr('r', 9.5)
            .attr('fill', '#000')
        } else if (d.shape === 'polygon') {
          legend.append('polygon')
            .attr('points', `${-19}, ${19} ${0}, ${19} ${-9.5}, ${0}`)
            .attr('fill', '#000')
        } else if (d.shape === 'rect') {
          legend.append('rect')
            .attr('x', -19)
            .attr('y', 0)
            .attr('width', 19)
            .attr('height', 19)
            .attr('fill', '#000')
        }
      })
)}

function _26(legendPfam,pfamShapes){return(
legendPfam.selectAll("text")
  	.data(pfamShapes)
  	.enter().append("text")
    	.attr("x", 0)
    	.attr("y", (d, i) => i * 20 + 9.5)
    	.attr("dy", "0.35em")
    	.text(d => d.pfam)
)}

function _legendHit(svg){return(
svg.append("g")
  .attr("transform", "translate(30, 400)")
)}

function _hitValues(){return(
[1, 5, 10]
)}

function _29(legendHit,hitValues,radius){return(
legendHit.selectAll('circle')
    .data(hitValues)
    .join('circle')
      .attr('cx', -9.5)
      .attr('cy', (d, i) => i * 20 + 6.5)
      .attr('r', d => radius(d))
      .attr('fill', 'none')
      .attr('stroke', '#000')
)}

function _30(legendHit,hitValues){return(
legendHit.selectAll("text")
  	.data(hitValues)
  	.enter().append("text")
    	.attr("x", 10)
    	.attr("y", (d, i) => i * 20 + 6.5)
    	.attr("dy", "0.35em")
    	.text(d => d)
)}

export default function define(runtime, observer) {
  const main = runtime.module();
  function toString() { return this.url; }
  const fileAttachments = new Map([
    ["countries-50m.json", {url: new URL("./files/f4afb2d49f0b38843f6d74521b33d41f371246e1acd674ed78016dca816cb1d262b7c54b95d395a4dad7fba5d58ed19db2944698360d19483399c79565806794.json", import.meta.url), mimeType: "application/json", toString}],
    ["map_addice@1.csv", {url: new URL("./files/280dd620445720bdf9ae73047dc37b59b01216c9fd30957ff983c417162bfacbeb78542107209ffa319a594a69cd6fbfd3ed6898aea3ebb2d856636f24ebfbf1.csv", import.meta.url), mimeType: "text/csv", toString}]
  ]);
  main.builtin("FileAttachment", runtime.fileAttachments(name => fileAttachments.get(name)));
  main.variable(observer()).define(["md"], _1);
  main.variable(observer("map_orign")).define("map_orign", ["htl","width","height","path","outline","graticule","land","borders"], _map_orign);
  main.variable(observer("borders")).define("borders", ["topojson","world"], _borders);
  main.variable(observer("projection")).define("projection", ["d3"], _projection);
  main.variable(observer("path")).define("path", ["d3","projection"], _path);
  main.variable(observer("outline")).define("outline", _outline);
  main.variable(observer("height")).define("height", ["d3","projection","width","outline"], _height);
  main.variable(observer("graticule")).define("graticule", ["d3"], _graticule);
  main.variable(observer("land")).define("land", ["topojson","world"], _land);
  main.variable(observer("svg")).define("svg", ["d3","map_orign"], _svg);
  main.variable(observer("world")).define("world", ["FileAttachment"], _world);
  main.variable(observer()).define(["svg","topojson","world","path"], _12);
  main.variable(observer()).define(["svg","topojson","world","path"], _13);
  main.variable(observer("map")).define("map", ["__query","FileAttachment","invalidation"], _map);
  main.variable(observer("radius")).define("radius", ["d3","map"], _radius);
  main.variable(observer("color")).define("color", ["d3","map"], _color);
  main.variable(observer("color_1")).define("color_1", ["d3"], _color_1);
  main.variable(observer("bubbles")).define("bubbles", ["svg","map","projection"], _bubbles);
  main.variable(observer()).define(["bubbles","d3","radius","color_1"], _19);
  main.variable(observer("legendHabitat")).define("legendHabitat", ["svg"], _legendHabitat);
  main.variable(observer()).define(["legendHabitat","color_1"], _21);
  main.variable(observer()).define(["legendHabitat","color_1"], _22);
  main.variable(observer("legendPfam")).define("legendPfam", ["svg"], _legendPfam);
  main.variable(observer("pfamShapes")).define("pfamShapes", _pfamShapes);
  main.variable(observer()).define(["legendPfam","pfamShapes","d3"], _25);
  main.variable(observer()).define(["legendPfam","pfamShapes"], _26);
  main.variable(observer("legendHit")).define("legendHit", ["svg"], _legendHit);
  main.variable(observer("hitValues")).define("hitValues", _hitValues);
  main.variable(observer()).define(["legendHit","hitValues","radius"], _29);
  main.variable(observer()).define(["legendHit","hitValues"], _30);
  return main;
}

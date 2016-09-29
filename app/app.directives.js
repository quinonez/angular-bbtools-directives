'use strict';

/* Directives */
var BBTDirectives = angular.module('BBTDirectives', []);


BBTDirectives.directive('ejemploD3', function() {
	var link = function(scope, element, attr){
	var color = d3.scale.category10();
	var data = scope.data;
	var width = 500;
	var height = 700;
	var min = Math.min(width, height);
	var p= d3.select(element[0]).append("p");
	p.text("Conjunto A");
	var svg = d3.select(element[0]).append("svg");
	svg.attr({width: width, height: height});
	var conjuntoDatos = d3.range(15).map(function() { return {radius: Math.random() * 20 + 4}; });
	var  fuerza = d3.layout.force()
    			.gravity(0.05)
    			.size([500, 150])
			.nodes(conjuntoDatos);
 svg.selectAll("circle")
			.data(conjuntoDatos)
			.enter()
			.append("circle")
			.attr("r", function(d){return d.radius;})
			.style("fill", function(d, i) {
			return color(i);
			})
      			.on("dblclick", dblclick)
			.call(fuerza.drag);

var root = conjuntoDatos[0];
root.fixed = true;
var drag = fuerza.drag()
    	.on("dragstart", dragstart);
fuerza.start();
fuerza.on("tick", function(e) {
  var q = d3.geom.quadtree(conjuntoDatos),
      i = 0,
      n = conjuntoDatos.length;

  while (++i < n) q.visit(collide(conjuntoDatos[i]));
	svg.selectAll("circle")
      		.attr("cx", function(d,i) { return d.x; })
      		.attr("cy", function(d,i) { return d.y; });

});
/*+svg.on("click", function() {
  var p1 = d3.mouse(this);
  console.log("arriba");
  console.log(p1[1]);
  fuerza.resume();
});*/
function collide(node) {
  var r = node.radius + 16,
      nx1 = node.x - r,
      nx2 = node.x + r,
      ny1 = node.y - r,
      ny2 = node.y + r;
  return function(quad, x1, y1, x2, y2) {
    if (quad.point && (quad.point !== node)) {
      var x = node.x - quad.point.x,
          y = node.y - quad.point.y,
          l = Math.sqrt(x * x + y * y),
          r = node.radius + quad.point.radius;
      if (l < r) {
        l = (l - r) / l * .5;
        node.x -= x *= l;
        node.y -= y *= l;
        quad.point.x += x;
        quad.point.y += y;
      }
    }
    return x1 > nx2 || x2 < nx1 || y1 > ny2 || y2 < ny1;
  };
}
function dragstart(d) {
  var pos = d3.select(this).classed("fixed", d.fixed = true);
//pos.attr("cx", function(d) { return this.cx; });
console.log(this.cx);
}
function dblclick(d) {
  d3.select(this).classed("fixed", d.fixed = false);
}
	//svg.attr({width: width, height: height});
	/*var g = svg.append('g')
	// center the donut chart
	.attr('transform', 'translate(' + width / 2 + ',' + height / 2 + ')');
	// add the <path>s for each arc slice
	g.selectAll('path').data(pie(data))
	.enter().append('path')
	.style('stroke', 'white')
	.attr('d', arc)
	.attr('fill', function(d, i){ return color(i) });*/
	/*var svg = d3.select(element[0]).selectAll("div")
	.data(data) // <-- The answer is here!
	.enter().append("svg");
	svg.attr({width: width, height: height});*/

	};

	return {
	link : link,
	restrict:'E',
	scope: { data: '=' }
	};    
	});





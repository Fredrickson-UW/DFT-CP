
var scene = new THREE.Scene();
var boxheight = document.getElementById("canvasbox").offsetHeight;
var boxwidth = document.getElementById("canvasbox").offsetWidth;
var aspectratio = boxwidth / boxheight;
//var camera = new THREE.PerspectiveCamera( 45, aspectratio, 1, 1000 );
var camera = new THREE.OrthographicCamera(-10*aspectratio,10*aspectratio,-10,10 , 10, 1000 );
var renderer = new THREE.WebGLRenderer({preserveDrawingBuffer: true});
var format = '.jpg'

renderer.setSize(boxwidth, boxheight, false);
renderer.domElement.id = "myvanvas";
window.addEventListener('resize', resizefunction);
function resizefunction() {
  boxheight = document.getElementById("canvasbox").offsetHeight;
  boxwidth = document.getElementById("canvasbox").offsetWidth;
  renderer.setSize( boxwidth, boxheight );
  aspectratio = boxwidth / boxheight;
  camera.left = -10*aspectratio;
  camera.right = 10*aspectratio;
  camera.updateProjectionMatrix();
}
document.getElementById("canvasbox").appendChild( renderer.domElement );
renderer.gammaInput = true;
renderer.gammaOutput = true;

// Jonathan has replaced the use of OrbitControls with TrackballControls, enabling rotation along all three axes
//var controls = new THREE.OrbitControls( camera, renderer.domElement);
//controls.update();
var controls = new THREE.TrackballControls( camera, renderer.domElement );
controls.rotateSpeed = 4;
controls.update();

var pos = new THREE.Vector3(0,2,3);
var light = new THREE.PointLight(0xffffff,0.5);
var amblight = new THREE.AmbientLight(0xffffff,0.5)
light.position.set(.1,.1,-30);
scene.add(light);
scene.add(amblight);
camera.position.set(0, 0, 25);
scene.add(camera);
scene.background = new THREE.Color(0xffffff);
const lighthelper = new THREE.PointLightHelper(light);
scene.add(lighthelper);

function animate() {
  requestAnimationFrame( animate );
  renderer.render( scene, camera );
  controls.update();
  lighthelper.update();
}
animate();

function define_material(color_rgb,alpha) {
  var new_material = new THREE.MeshPhysicalMaterial(color_rgb);
	if(alpha < 1) {
    new_material.transparent = true;
	}
	new_material.side = THREE.DoubleSide;
	new_material.opacity = alpha;
//new_material.envMap = THREE.SphericalReflectionMapping;
  new_material.reflectivity = 10;
  new_material.roughness = 0;
  new_material.metalness = 0;
//new_material.specular = new THREE.Color(0xffffff);
  return(new_material);
}

function define_material2(color_rgb,alpha) {
  var new_material = new THREE.MeshPhysicalMaterial(color_rgb);
	if(alpha < 1) {
    new_material.transparent = true;
	}
	new_material.side = THREE.BackSide;
	new_material.opacity = alpha;
//new_material.envMap = THREE.SphericalReflectionMapping;
  new_material.reflectivity = 10;
  new_material.roughness = 0;
  new_material.metalness = 0;
//new_material.specular = new THREE.Color(0xffffff);
  return(new_material);
}

function add_sphere(x,y,z,material,radius) {
  var geometry2 = new THREE.SphereGeometry(radius,40,40);
  geometry2.translate(x,y,z)
  geometry2.computeVertexNormals();
  var new_sphere2 = new THREE.Mesh( geometry2, material);
  scene.add( new_sphere2) ;
  return(new_sphere2);
}

function drawCPs(xyz,abc,template,CPs,scale) {
  var natoms_xyz = xyz.length;
  var natoms_template = template.length;
  var lobes = [];
  var lobes_temp = [];
  var counter = 0;
  var j1,j2,j3,k1,k2,x1,y1,z1,x2,y2,z2,dist;
  for( j1=-4; j1 < 5; j1++) {
    for( j2=-4; j2 < 5; j2++) {
      for( j3=-4; j3 < 5; j3++) {
        for(k1 = 0; k1 < natoms_xyz; k1++) {
          for(k2 = 0; k2 < natoms_template; k2++) {
            x1=xyz[k1][1]+abc[0][0]*j1+abc[1][0]*j2+abc[2][0]*j3;
            y1=xyz[k1][2]+abc[0][1]*j1+abc[1][1]*j2+abc[2][1]*j3;
            z1=xyz[k1][3]+abc[0][2]*j1+abc[1][2]*j2+abc[2][2]*j3;
            x2=template[k2][1];
            y2=template[k2][2];
            z2=template[k2][3];
            dist = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**0.5;
            if(dist < 0.1) {
              lobes_temp =  add_lobesCP(x2,y2,z2,plus,minus,CPs[k1],scale);
              lobes[counter]=lobes_temp[0];
              counter++;
              lobes[counter]=lobes_temp[1];
              counter++;
            }
          }
        }
      }
    }
  }
  return(lobes);
}


function add_lobes(x_o,y_o,z_o,bluematerial,greenmaterial,coeff,scale) {
  var lobes = [];
  var r = 0;
  var x,y,z,theta,phi;
  positive_surf = new THREE.SphereGeometry(1,60,60);
  negative_surf = new THREE.SphereGeometry(1,60,60);
  for( let j = 0; j < positive_surf.vertices.length; j++) {
    r=0;
    x=positive_surf.vertices[j].x;
	  y=positive_surf.vertices[j].y;
	  z=positive_surf.vertices[j].z;
    theta=Math.acos(z);
	  phi =0.0;
	  if((theta > 0) && (theta < Math.PI)) {
      phi = Math.acos(x/Math.sin(theta));
      if(x/Math.sin(theta) > 1) {
        phi = 0;
      }
      if(x/Math.sin(theta) < -1) {
        phi = Math.PI;
      }
    }
    if(y < 0.0) {
      phi=-phi;
    }
    r += scale*0.5*(1/Math.PI)**0.5*coeff[0];
    r += scale*(3/(4*Math.PI))**0.5*Math.sin(theta)*Math.cos(phi)*coeff[1];
    r += scale*(3/(4*Math.PI))**0.5*Math.sin(theta)*Math.sin(phi)*coeff[2];
	  r += scale*(3/(4*Math.PI))**0.5*Math.cos(theta)*coeff[3];
	  r += scale*0.25*(15/(Math.PI))**0.5*Math.sin(theta)**2*(Math.cos(phi)**2-Math.sin(phi)**2)*coeff[4];
	  r += scale*0.25*(5/(Math.PI))**0.5*(-(Math.sin(theta)**2)*(Math.cos(phi)**2+Math.sin(phi)**2)+2*(Math.cos(theta))**2)*coeff[5];
	  r += scale*0.5*(15/(Math.PI))**0.5*Math.sin(theta)**2*(Math.cos(phi)*Math.sin(phi))*coeff[6];
    r += scale*0.5*(15/(Math.PI))**0.5*Math.sin(theta)*Math.cos(phi)*Math.cos(theta)*coeff[7];
	  r += scale*0.5*(15/(Math.PI))**0.5*Math.sin(theta)*Math.sin(phi)*Math.cos(theta)*coeff[8];
    if (r<0) {
      positive_surf.vertices[j].x=0.0;
		  positive_surf.vertices[j].y=0.0;
		  positive_surf.vertices[j].z=0.0;
      negative_surf.vertices[j].x=-r*Math.sin(theta)*Math.cos(phi);
		  negative_surf.vertices[j].y=-r*Math.sin(theta)*Math.sin(phi);
		  negative_surf.vertices[j].z=-r*Math.cos(theta);
    }
    else {
      negative_surf.vertices[j].x=0.0;
      negative_surf.vertices[j].y=0.0;
      negative_surf.vertices[j].z=0.0;
      positive_surf.vertices[j].x=r*Math.sin(theta)*Math.cos(phi);
      positive_surf.vertices[j].y=r*Math.sin(theta)*Math.sin(phi);
      positive_surf.vertices[j].z=r*Math.cos(theta);
    }
  }
  positive_surf.translate(x_o,y_o,z_o)
  negative_surf.translate(x_o,y_o,z_o)
  positive_surf.verticesNeedUpdate = true;
  negative_surf.verticesNeedUpdate = true;
  positive_surf.facesNeedUpdate = true;
  negative_surf.facesNeedUpdate = true;
  positive_surf.elementsNeedUpdate = true;
  negative_surf.uvsNeedUpdate = true;
  positive_surf.computeFaceNormals();
  negative_surf.computeFaceNormals();
  positive_surf.computeVertexNormals();
  negative_surf.computeVertexNormals();
  lobes[0] = new THREE.Mesh(positive_surf, bluematerial);
  lobes[1] = new THREE.Mesh(negative_surf, greenmaterial);
  scene.add(lobes[0]);
  scene.add(lobes[1]);
  return(lobes);
}

/*
if(voxel_r > 0.0) {
  for(l=1;l<lmax+1;l++) {
    m=0;
    XSFOUT->int_Ylm[central_atom][l][m]+=1*gsl_sf_legendre_sphPlm(l,m,costheta)*XSFIN->grid[jx][jy][jz]/vmap->neighbor_count[jx][jy][jz];
    for(m=1;m<l+1;m++) {
      XSFOUT->int_Ylm[central_atom][l][2*m-1]+=pow(2,0.5)*cosmphi[m]*gsl_sf_legendre_sphPlm(l,m,costheta)*XSFIN->grid[jx][jy][jz]/vmap->neighbor_count[jx][jy][jz];
      XSFOUT->int_Ylm[central_atom][l][2*m]+=pow(2,0.5)*sinmphi[m]*gsl_sf_legendre_sphPlm(l,m,costheta)*XSFIN->grid[jx][jy][jz]/vmap->neighbor_count[jx][jy][jz];
    }
  }
}
*/

function add_lobesCP(x_o,y_o,z_o,bluematerial,greenmaterial,coeff,scale) {
  var lobes = [];
  var r = 0;
  var x,y,z,theta,phi;
  positive_surf = new THREE.SphereGeometry(1,40,40);
  negative_surf = new THREE.SphereGeometry(1,40,40);
  for( let j = 0; j < positive_surf.vertices.length; j++) {
    r=0;
    x=positive_surf.vertices[j].x;
	  y=positive_surf.vertices[j].y;
	  z=positive_surf.vertices[j].z;
    theta=Math.acos(z);
	  phi =0.0;
	  if((theta > 0) && (theta < Math.PI)) {
      phi = Math.acos(x/Math.sin(theta));
		  if(x/Math.sin(theta) > 1) {
        phi = 0;
		  }
	    if(x/Math.sin(theta) < -1) {
        phi = Math.PI;
		  }
    }
    if(y < 0.0) {
      phi=-phi;
    }
    r += scale*0.5*(1/Math.PI)**0.5*coeff[0];
	  r += scale*0.5*(3/(Math.PI))**0.5*Math.cos(theta)*coeff[1];
	  r += -scale*(3/(8*Math.PI))**0.5*Math.sin(theta)*Math.cos(phi)*(2)**0.5*coeff[2];
    r += -scale*(3/(8*Math.PI))**0.5*Math.sin(theta)*Math.sin(phi)*(2)**0.5*coeff[3];
    r += scale*0.25*(5/(Math.PI))**0.5*(3*Math.cos(theta)**2-1)*coeff[4];
    r += -scale*0.5*(15/(Math.PI))**0.5*Math.sin(theta)*Math.cos(phi)*Math.cos(theta)*coeff[5];
	  r += -scale*0.5*(15/(Math.PI))**0.5*Math.sin(theta)*Math.sin(phi)*Math.cos(theta)*coeff[6];
	  r += scale*0.25*(15/(Math.PI))**0.5*Math.sin(theta)**2*Math.cos(2*phi)*coeff[7];
	  r += scale*0.25*(15/(Math.PI))**0.5*Math.sin(theta)**2*Math.sin(2*phi)*coeff[8];
    r += scale*0.25*(7/(Math.PI))**0.5*(5*Math.cos(theta)**3-3*Math.cos(theta))*coeff[9];
	  r += -scale*0.125*(21/(Math.PI))**0.5*Math.sin(theta)*(5*Math.cos(theta)**2-1.0)*Math.cos(phi)*(2)**0.5*coeff[10];
	  r += -scale*0.125*(21/(Math.PI))**0.5*Math.sin(theta)*(5*Math.cos(theta)**2-1.0)*Math.sin(phi)*(2)**0.5*coeff[11];
	  r += scale*0.25*(105/(2*Math.PI))**0.5*Math.sin(theta)**2*(Math.cos(theta))*Math.cos(2*phi)*(2)**0.5*coeff[12];
	  r += scale*0.25*(105/(2*Math.PI))**0.5*Math.sin(theta)**2*(Math.cos(theta))*Math.sin(2*phi)*(2)**0.5*coeff[13];
	  r += -scale*0.125*(35/(Math.PI))**0.5*Math.sin(theta)**3*Math.cos(3*phi)*(2)**0.5*coeff[14];
	  r += -scale*0.125*(35/(Math.PI))**0.5*Math.sin(theta)**3*Math.sin(3*phi)*(2)**0.5*coeff[15];

	  r += scale*(3/16)*(1/(Math.PI))**0.5*(35*Math.cos(theta)**4-30*Math.cos(theta)**2+3)*coeff[16];
	  r += -scale*(3/8)*(5/(Math.PI))**0.5*Math.sin(theta)*(7*Math.cos(theta)**3-3*Math.cos(theta))*Math.cos(phi)*(2)**0.5*coeff[17];
	  r += -scale*(3/8)*(5/(Math.PI))**0.5*Math.sin(theta)*(7*Math.cos(theta)**3-3*Math.cos(theta))*Math.sin(phi)*(2)**0.5*coeff[18];
	  r += scale*(3/8)*(5/(2*Math.PI))**0.5*Math.sin(theta)**2*(7*Math.cos(theta)**2-1)*Math.cos(2*phi)*(2)**0.5*coeff[19];
	  r += scale*(3/8)*(5/(2*Math.PI))**0.5*Math.sin(theta)**2*(7*Math.cos(theta)**2-1)*Math.sin(2*phi)*(2)**0.5*coeff[20];
	  r += -scale*(3/8)*(35/(Math.PI))**0.5*Math.sin(theta)**3*Math.cos(theta)*Math.cos(3*phi)*(2)**0.5*coeff[21];
	  r += -scale*(3/8)*(35/(Math.PI))**0.5*Math.sin(theta)**3*Math.cos(theta)*Math.sin(3*phi)*(2)**0.5*coeff[22];
	  r += scale*(3/16)*(35/(2*Math.PI))**0.5*Math.sin(theta)**4*Math.cos(4*phi)*(2)**0.5*coeff[23];
	  r += scale*(3/16)*(35/(2*Math.PI))**0.5*Math.sin(theta)**4*Math.sin(4*phi)*(2)**0.5*coeff[24];

    if (r<0) {
      positive_surf.vertices[j].x=0.0;
		  positive_surf.vertices[j].y=0.0;
	    positive_surf.vertices[j].z=0.0;
      negative_surf.vertices[j].x=-r*Math.sin(theta)*Math.cos(phi);
		  negative_surf.vertices[j].y=-r*Math.sin(theta)*Math.sin(phi);
		  negative_surf.vertices[j].z=-r*Math.cos(theta);
	  }
	  else {
      negative_surf.vertices[j].x=0.0;
	    negative_surf.vertices[j].y=0.0;
	    negative_surf.vertices[j].z=0.0;
	    positive_surf.vertices[j].x=r*Math.sin(theta)*Math.cos(phi);
	    positive_surf.vertices[j].y=r*Math.sin(theta)*Math.sin(phi);
	    positive_surf.vertices[j].z=r*Math.cos(theta);
	  }
  }
  positive_surf.translate(x_o,y_o,z_o)
  negative_surf.translate(x_o,y_o,z_o)
  positive_surf.verticesNeedUpdate = true;
  negative_surf.verticesNeedUpdate = true;
  positive_surf.facesNeedUpdate = true;
  negative_surf.facesNeedUpdate = true;
  positive_surf.elementsNeedUpdate = true;
  negative_surf.uvsNeedUpdate = true;
  positive_surf.computeFaceNormals();
  negative_surf.computeFaceNormals();
  positive_surf.computeVertexNormals();
  negative_surf.computeVertexNormals();
  lobes[0] = new THREE.Mesh(positive_surf, bluematerial);
  lobes[1] = new THREE.Mesh(negative_surf, greenmaterial);
  scene.add(lobes[0]);
  scene.add(lobes[1]);
  return(lobes);
}

function drawbonds_cylinder(atoms,atom1,color1,atom2,color2,d_min,d_max,radius) {
  var dist = 0.0;
  var d_x,d_y,d_z,theta,phi;
  var new_cylinder = [];
  var counter = 0;
  var nbonds=0;
  for(let j = 0; j < atoms.length ; j++) {
    if(atom1 === atom2) {
      for(let k = 0; k < atoms.length ; k++) {
        dist = ((atoms[j][1]-atoms[k][1])**2.0 + (atoms[j][2]-atoms[k][2])**2.0 + (atoms[j][3]-atoms[k][3])**2.0)**0.5;
        if((dist >= d_min) && (dist <= d_max) && (atoms[j][0] === atom1) && (atoms[k][0] === atom2)) {
          nbonds++;
	        d_x = (atoms[k][1]-atoms[j][1]);
	        d_y = (atoms[k][2]-atoms[j][2]);
	        d_z = (atoms[k][3]-atoms[j][3]);
	        theta=Math.acos(d_z/dist);
	        phi =0.0;
	        if((theta > 0) && (theta < Math.PI)) {
            var acin = Math.round( d_x / Math.sin(theta) / dist * 1000000 + Number.EPSILON ) / 1000000;
            phi = Math.acos(acin);
	        }
	        if(d_y < 0.0) {
	          phi=-phi;
	        }

	        var geometry2 = new THREE.CylinderGeometry(radius,radius,dist,20);
	        geometry2.rotateX(Math.PI/2);
	        geometry2.rotateY(theta);
	        geometry2.rotateZ(phi);
	        geometry2.translate(atoms[j][1]+d_x*0.5,atoms[j][2]+d_y*0.5,atoms[j][3]+d_z*0.5);
	        geometry2.computeFaceNormals();
	        geometry2.normalsNeedUpate = true;
	        new_cylinder[counter] = new THREE.Mesh( geometry2,color1);
	        scene.add(new_cylinder[counter]);
	        counter++;
        }
      }
    }
    else {
      for(let k = 0; k < atoms.length ; k++) {
        dist = ((atoms[j][1]-atoms[k][1])**2.0 + (atoms[j][2]-atoms[k][2])**2.0 + (atoms[j][3]-atoms[k][3])**2.0)**0.5;
        if((dist >= d_min) && (dist <= d_max) && (atoms[j][0] === atom1) && (atoms[k][0] === atom2)) {
          d_x = (atoms[k][1]-atoms[j][1]);
	        d_y = (atoms[k][2]-atoms[j][2]);
		      d_z = (atoms[k][3]-atoms[j][3]);
		      theta=Math.acos(d_z/dist);
		      phi =0.0;
		      if((theta > 0) && (theta < Math.PI)) {
            phi = Math.acos(d_x/Math.sin(theta)/dist);
		      }
		      if(d_y < 0.0) {
            phi=-phi;
		      }
	        var geometry2 = new THREE.CylinderGeometry(radius,radius,dist/2.0,40);
          geometry2.rotateX(Math.PI/2);
	        geometry2.rotateY(theta);
          geometry2.rotateZ(phi);
		      geometry2.translate(atoms[j][1]+d_x*0.25,atoms[j][2]+d_y*0.25,atoms[j][3]+d_z*0.25);
		      geometry2.computeFaceNormals();
		      geometry2.normalsNeedUpate = true;
		      new_cylinder[counter] = new THREE.Mesh( geometry2,color1);
		      scene.add(new_cylinder[counter]);
		      counter++;
       	  geometry2 = new THREE.CylinderGeometry(radius,radius,dist/2.0,40);
		      geometry2.rotateX(Math.PI/2);
		      geometry2.rotateY(theta);
		      geometry2.rotateZ(phi);
		      geometry2.translate(atoms[j][1]+d_x*0.75,atoms[j][2]+d_y*0.75,atoms[j][3]+d_z*0.75);
		      geometry2.computeFaceNormals();
		      geometry2.normalsNeedUpate = true;
		      new_cylinder[counter] = new THREE.Mesh( geometry2,color2);
		      scene.add(new_cylinder[counter]);
		      counter++;
        }
      }
    }
  }
  return(new_cylinder);
}

function center_controls(atoms) {
  var x,y,z;
  x=0;
  y=0;
  z=0;
  for (let j = 0; j < atoms.length ; j++) {
    x+=atoms[j][1]/atoms.length;
    y+=atoms[j][2]/atoms.length;
    z+=atoms[j][3]/atoms.length;
  }
  controls.target = new THREE.Vector3(x,y,z);
}

function drawmol_spheres(atoms,atom,material,radius) {
  var new_sphere = [];
  var counter=0;
  for (let j = 0; j < atoms.length ; j++) {
    if(atoms[j][0]===atom) {
      new_sphere[counter] = add_sphere(atoms[j][1],atoms[j][2],atoms[j][3],material,radius);
	    counter++
    }
  }
  return(new_sphere);
}

function draw_plane(material) {
  var geometry = new THREE.Geometry();
  geometry.vertices.push(
  	new THREE.Vector3( -1, -1,  0 ),
  	new THREE.Vector3( -1,  1,  0 ),
  	new THREE.Vector3(  1, -1,  0 )
  );
  geometry.faces.push(new THREE.Face3( 0, 1, 2 ) );
  geometry.computeBoundingSphere();
  var mesh = new THREE.Mesh( geometry, material );
  scene.add( mesh );
}

function draw_polygon(geometry, material) {
  var mesh = new THREE.Mesh( geometry, material );
  scene.add( mesh);
}

function remove_objects(spheres) {
  for (let j = 0; j < spheres.length ; j++) {
    scene.remove(spheres[j]);
  }
}

function add_objects(spheres) {
  for (let j = 0; j < spheres.length ; j++) {
    scene.add(spheres[j]);
  }
}

function dontshow_objects(spheres) {
  for (let j=0; j < spheres.length; j++) {
    spheres[j].visible = false;
  }
}

function show_objects(spheres) {
  for (let j=0; j < spheres.length; j++) {
    spheres[j].visible = true;
  }
}

function hide_objects(spheres) {
  for (let j = 0; j < spheres.length ; j++) {
    if(spheres[j].visible === false) {
      spheres[j].visible = true;
    }
    else {
      spheres[j].visible = false;
    }
  }
}

function hextorgb(hex) {
  var shorthandRegex = /^#?([a-f\d])([a-f\d])([a-f\d])$/i;
  hex = hex.replace(shorthandRegex, function(m, r, g, b) {
    return r + r + g + g + b + b;
  });
  var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  return result ? {
    r: parseInt(result[1], 16),
    g: parseInt(result[2], 16),
    b: parseInt(result[3], 16)
  } : null;
}

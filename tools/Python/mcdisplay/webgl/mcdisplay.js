var transformPoints = function(apoints, transform) {
    // Create a BufferGeometry
    var geometry = new THREE.BufferGeometry();

    // Convert apoints array to a Float32Array and set it as the position attribute
    var vertices = new Float32Array(apoints.length * 3);
    for (var i = 0; i < apoints.length; i++) {
        vertices[i * 3] = apoints[i].x;
        vertices[i * 3 + 1] = apoints[i].y;
        vertices[i * 3 + 2] = apoints[i].z;
    }
    geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));

    // Apply the transformation matrix
    geometry.applyMatrix4(transform);

    // Extract the transformed vertices
    var transformedVertices = [];
    var position = geometry.getAttribute('position');
    for (var i = 0; i < position.count; i++) {
        transformedVertices.push(new THREE.Vector3(position.getX(i), position.getY(i), position.getZ(i)));
    }

    return transformedVertices;
}

// "main" class which is a collection of scene graph setup routives and data objects
//
var Main = function ()
{
    this.scene;
    this.camera;
    this.renderer;
    this.controls;

    this.compnodes = {};
    this.rootnode = new THREE.Object3D();
    this.bbnode = new THREE.Object3D();
    this.rootnode.add(this.bbnode);

    this.compColors = [0xf80000, 0x00f800, 0x0000f8, 0x00f8f8, 0xf800f8, 0x00f880, 0xf8f800, 0xf88000, 0x80f800, 0x0080f8, 0x800f8, 0xf80080, 0xa8a8a8];
    this.iColor = -1;

    this.rayColor = 0x00ffff;
    this.aCompVertices
    this.rayobj;
    this.raynodes = [];

    this.iRay = -1;

    this.lut = null;
}
// add circle to xy, xz or yz plane
//
Main.prototype.addCircle = function(plane, x, y, z, radius, parent, linecolor)
{
    if (radius == 0)
    {
        return;
    }

    if (plane === 'xy') {
        Main.prototype.addNewCircle(x, y, z, radius, 0, 0, 1, 32, parent, linecolor)
    }
    if (plane === 'xz') {
        Main.prototype.addNewCircle(x, y, z, radius, 0, 1, 0, 32, parent, linecolor)
    }
    if (plane === 'yz') {
        Main.prototype.addNewCircle(x, y, z, radius, 1, 0, 0, 32, parent, linecolor)
    }
}
// add sphere
//		x,y,z - center coordinates
//		radius
//		wseg 	- width segments
//		hseg 	- height segments
Main.prototype.addSphere = function(x, y, z, radius, wseg, hseg, parent, color)
{
    if (radius === 0)
    {
        return;
    }
    var geometry = new THREE.SphereGeometry(radius, wseg, hseg);
    var material = new THREE.MeshLambertMaterial( {color: color} );
    var sphere = new THREE.Mesh( geometry, material );

    sphere.position.x = x;
    sphere.position.y = y;
    sphere.position.z = z;

    parent.add( sphere );
}

Main.prototype.addCone = function(x, y, z, radius, height, nx, ny, nz, radSeg, parent, color)
{
    if (radius === 0)
    {
        return;
    }
    var geometry = new THREE.ConeGeometry(radius, height, radSeg);
    var material = new THREE.MeshLambertMaterial( {color: color} );
    var cone = new THREE.Mesh( geometry, material );

    cone.position.x = x;
    cone.position.y = y;
    cone.position.z = z;

    var original_axis = new THREE.Vector3(0, 1, 0);
    var align_axis = new THREE.Vector3(nx, ny, nz);
    cone.quaternion.setFromUnitVectors(original_axis, align_axis.clone().normalize());

    parent.add( cone );
}

Main.prototype.addCylinder = function(x, y, z, radius, height, thickness, nx, ny, nz, radSeg, parent, color)
{
    let geometry = new THREE.CylinderGeometry(radius, radius, height, radSeg);
    let original_axis = new THREE.Vector3(0, 1, 0);

    if(thickness > 0){
        //correcting some differences between ExtrudeGeometry and CylinderGeometry implementation
        original_axis = new THREE.Vector3(0, 0, 1);
        y = -height/2;

        let outer_cylinder = new THREE.Shape();
        outer_cylinder.absarc(0, 0, radius, 0, Math.PI * 2);
        let inner_cylinder = new THREE.Path();
        inner_cylinder.absarc(0, 0, (radius - thickness), 0, Math.PI * 2);
        outer_cylinder.holes.push(inner_cylinder);

        geometry = new THREE.ExtrudeGeometry(outer_cylinder, {
            steps: 1,
            depth: height,
            curveSegments: radSeg,
            bevelEnabled: false
        });
    }


    const material = new THREE.MeshLambertMaterial({color: color});
    const cylinder = new THREE.Mesh(geometry, material);

    cylinder.position.x = x;
    cylinder.position.y = y;
    cylinder.position.z = z;

    var align_axis = new THREE.Vector3(nx, ny, nz);
    cylinder.quaternion.setFromUnitVectors(original_axis, align_axis.clone().normalize());

    parent.add( cylinder );
}

Main.prototype.addAnnulus = function(x, y, z, outer_radius, inner_radius, nx, ny, nz, radSeg, parent, color)
{
    let geometry = new THREE.RingGeometry(outer_radius-inner_radius, outer_radius, radSeg);
    let original_axis = new THREE.Vector3(0, 0, 1);

    const material = new THREE.MeshLambertMaterial({color: color, side: THREE.DoubleSide});
    const annulus = new THREE.Mesh(geometry, material);

    annulus.position.x = x;
    annulus.position.y = y;
    annulus.position.z = z;

    var align_axis = new THREE.Vector3(nx, ny, nz);
    annulus.quaternion.setFromUnitVectors(original_axis, align_axis.clone().normalize());

    parent.add( annulus );
}

Main.prototype.addDisc = function(x, y, z, radius, nx, ny, nz, radSeg, parent, color)
{
    Main.prototype.addAnnulus(x, y, z, radius, radius, nx, ny, nz, radSeg, parent, color);
}

Main.prototype.addNewCircle = function(x, y, z, radius, nx, ny, nz, radSeg, parent, color)
{
    Main.prototype.addAnnulus(x, y, z, radius, 0.01, nx, ny, nz, radSeg, parent, color);
}

Main.prototype.addBox = function(x, y, z, a, b, c, thickness, parent, color)
{
    let geometry = new THREE.BoxGeometry(a, b, c);

    if(thickness > 0){
        a /= 2;
        b /= 2;

        let outerBox = new THREE.Shape();
        outerBox.moveTo(-a, -b); // Bottom-left corner
        outerBox.lineTo(a, -b);  // Bottom-right corner
        outerBox.lineTo(a, b);   // Top-right corner
        outerBox.lineTo(-a, b);  // Top-left corner

        const a_inner = a - thickness;
        const b_inner = b - thickness;

        let innerBox =  new THREE.Shape();
        innerBox.moveTo(-a_inner, -b_inner);
        innerBox.lineTo(a_inner, -b_inner);
        innerBox.lineTo(a_inner, b_inner);
        innerBox.lineTo(-a_inner, b_inner);

        outerBox.holes.push(innerBox);

        geometry = new THREE.ExtrudeGeometry(outerBox, {
            steps: 1,
            depth: c,
            bevelEnabled: false
        });
    }

    const material = new THREE.MeshLambertMaterial({color: color});
    const box = new THREE.Mesh(geometry, material);

    box.position.x = x;
    box.position.y = y;
    box.position.z = z;

    parent.add( box );
}

Main.prototype.addPolygon = function(faces, vertices, parent, color)
{
    let geometry = new THREE.BufferGeometry();

    //check faces for rank 4 TODO
    if(true){
        //transform faces to rank 3
        faces = getRank3Indices(faces);
    }

    geometry.setIndex(faces);
    geometry.setAttribute( 'position', new THREE.BufferAttribute( vertices, 3 ) );

    const material = new THREE.MeshLambertMaterial({color: color});
    const polygon = new THREE.Mesh(geometry, material);

    parent.add( polygon );
}

function getRank3Indices(rank4indices)
{
    let rank3indices = [];

    for(let i = 0; i < rank4indices.length; i+=4){
        let face = [];
        for(let j = 0; j < 4; j++){
            face.push(rank4indices[i+j]);
        }
        rank3indices.push(face[0]);
        rank3indices.push(face[1]);
        rank3indices.push(face[2]);

        rank3indices.push(face[0]);
        rank3indices.push(face[2]);
        rank3indices.push(face[3]);
    }
    return rank3indices;
}


// add pointlight
// 		center 	- THREE.Vector3 instance
Main.prototype.addLight = function(center)
{
    // pointlight (required to light up sphereMaterial)
    var pointLight = new THREE.PointLight(0xFFFFFF);
    pointLight.position.x = center.x;
    pointLight.position.y = center.y;
    pointLight.position.z = center.z;
    this.scene.add(pointLight);
}
//  initialize the scene
//
Main.prototype.init = function(campos, invert)
{
    this.scene = new THREE.Scene();

    this.camera = new THREE.PerspectiveCamera(60, window.innerWidth/window.innerHeight, 0.001, 130);
    this.camera.position.x = campos.x; // -50;
    this.camera.position.y = campos.y; // 0;
    this.camera.position.z = campos.z; // 50;

    this.renderer = new THREE.WebGLRenderer({
        antialias: true
    });
    this.renderer.setSize(window.innerWidth, window.innerHeight);
    if (invert)
    {
        this.renderer.setClearColor( 0xffffff );
    }

    document.body.appendChild(this.renderer.domElement);
    window.addEventListener("resize", event => {
        this.camera.aspect = window.innerWidth / window.innerHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(window.innerWidth, window.innerHeight);
    })

    // NOTE: initial camera view direction is along the x axis

    element = document.getElementById("3dcanvas");
    console.log(element);
    document.getElementById("3dcanvas").appendChild(this.renderer.domElement);

    this.controls = new THREE.OrbitControls(this.camera);
    this.controls.target.x = -campos.x/2; //1;
    this.controls.target.y = 0; //0;
    this.controls.target.z = campos.z; //49;

    this.camera.lookAt(this.controls.target);

    let light = new THREE.DirectionalLight(0xffffff, 0.5);
    light.position.setScalar(1);
    this.scene.add(light, new THREE.AmbientLight(0xffffff, 0.5));

    var gridXZ = new THREE.GridHelper(100, 100);
    this.scene.add(gridXZ);

    var gridXY = new THREE.GridHelper(100, 100);
    gridXY.rotation.x = Math.PI / 2;
    this.scene.add(gridXY);

    var gridYZ = new THREE.GridHelper(100, 100);
    gridYZ.rotation.z = Math.PI / 2;
    this.scene.add(gridYZ);


    this.scene.add(this.rootnode);
}
//  cat camera view according to campos, a Vector3
//
Main.prototype.setCameraView = function(campos)
{
    this.camera.position.x = campos.x;
    this.camera.position.y = campos.y;
    this.camera.position.z = campos.z;

    this.controls.target.x = 0;
    this.controls.target.y = 0;
    this.controls.target.z = campos.z;

    this.camera.lookAt(this.controls.target);
}

// add multiline
// 		arrVector3  -  an array of THREE.Vector3 instances
Main.prototype.addMultiLineV3 = function(arrVector3, parent, linecolor)
{
    var multilinematerial = new THREE.LineBasicMaterial({color: linecolor});
    var multilinegeometry = new THREE.BufferGeometry().setFromPoints(arrVector3);
    var multiline = new THREE.Line(multilinegeometry, multilinematerial);
    parent.add(multiline);
}
// add multiline proxy method
//		points  -  array of single points
Main.prototype.addMultiLine = function(points, parent, linecolor)
{
    vectors = [];
    for (var i = 0; i < points.length / 3; i++)
    {
        const v = new THREE.Vector3(points[i*3], points[i*3+1], points[i*3+2]);
        vectors.push(v);
    }
    this.addMultiLineV3(vectors, parent, linecolor);
}

//  adds tiny boxes to each particle ray scatter point
//
Main.prototype.putScatterPoints = function(raynode)
{
    var v;
    var vtx = raynode.vtx;
    for (var i = 0; i < vtx.length; i++)
    {
        var geometry = new THREE.BoxGeometry( 0.007, 0.007, 0.007 );
        v = vtx[i];
        geometry.translate(v.x, v.y, v.z);
        var material = new THREE.MeshBasicMaterial({ color: raynode.lut_color });
        var cube = new THREE.Mesh( geometry, material );

        raynode.annotcubes.add( cube );
    }
}

//  set color lookuptable min/max range
//  color maps: rainbow, cooltowarm, blackbody, grayscale
Main.prototype.setLutRange = function(min, max)
{
    this.lut = new THREE.Lut( "cooltowarm", 512 );
    this.lut.setMin(min);
    this.lut.setMax(max);
}
// add a ray node
//
Main.prototype.addRayNode = function(rayObj, vertices, speed)
{
    var lut_color = this.lut.getColor(speed);
    rayObj.lut_color = lut_color;

    var multilinematerial = new THREE.LineBasicMaterial({color: lut_color});
    var multilinegeometry = new THREE.BufferGeometry().setFromPoints(vertices);
    var multiline = new THREE.Line(multilinegeometry, multilinematerial);
    rayObj.add(multiline);

    rayObj.vtx = vertices; // WARNING: this is a hidden field hacky stuff, but only used in the function putScatterPoints
    annotcubes = new THREE.Object3D();
    rayObj.add(annotcubes);
    rayObj.annotcubes = annotcubes;

    this.raynodes.push(rayObj);
    rayObj.visible = false;
}
// add a component origo
//		compname  -  component name for the reference dict
Main.prototype.addComponent = function(compname)
{
    var comp = new THREE.Object3D();
    this.rootnode.add(comp);
    this.compnodes[compname] = comp;
    return comp;
}
// returns a cyclically new color intended for a component draw node
//
Main.prototype.getNextComponentColor = function()
{
    this.iColor += 1;
    if (this.iColor >= this.compColors.length)
    {
        this.iColor = 0;
    }
    return this.compColors[this.iColor];
}
// show all rays in this.raynodes
//
Main.prototype.showAllRays = function()
{
    for (var i = 0; i < this.raynodes.length; i++)
    {
        node = this.raynodes[i];
        // NOTE: the following two lines should instead hvae been implemented via a function called "update" and set in a proper manor
        if (node.parent == null) { this.putScatterPoints(node); this.rootnode.add(node); };
        this.raynodes[i].visible = true;
    }
}
// remove all rays in this.raynodes from the scene graph
//
Main.prototype.hideAllRays = function()
{
    for (var i = 0; i < this.raynodes.length; i++)
    {
        this.raynodes[i].visible = false;
    }
}
// hide all annotation cubes on all rays
//
Main.prototype.hideCubes = function()
{
    for (var i = 0; i < this.raynodes.length; i++)
    {
        ray = this.raynodes[i];
        ray.annotcubes.visible = false;
    }
}
// show all annotation cubes on all rays
//
Main.prototype.showCubes = function()
{
    for (var i = 0; i < this.raynodes.length; i++)
    {
        ray = this.raynodes[i];
        ray.annotcubes.visible = true;
    }
}
// iterates to attach the next ray in the global ray sequence
//
Main.prototype.showNextRay = function()
{
    if (this.raynodes.length == 0) { return; }

    lastRay = this.iRay;
    this.iRay += 1;
    if (this.iRay >= this.raynodes.length)
    {
        this.iRay = 0;
    }
    newnode = this.raynodes[this.iRay];
    // NOTE: the following two lines should instead hvae been implemented via a function called "update" and set in a proper manor
    if (node.parent == null) { this.putScatterPoints(newnode); this.rootnode.add(newnode); };
    newnode.visible = true;

    this.raynodes[lastRay].visible = false;
}
// iterates to attach the next ray in the global ray sequence
//
Main.prototype.hideRay = function(idx)
{
    if (this.raynodes.length == 0) { return idx; }

    if (idx >= this.raynodes.length | idx < 0)
    {
        throw "idx out of range"
    }
    this.raynodes[idx].visible = false;
}
// iterates to attach the next ray in the global ray sequence
//
Main.prototype.showRay = function(idx)
{
    if (this.raynodes.length == 0) { return idx; }

    if (idx >= this.raynodes.length | idx < 0)
    {
        throw "idx out of range"
    }
    //this.rootnode.add(this.raynodes[idx]);
    node = this.raynodes[idx];
    // NOTE: the following two lines should instead hvae been implemented via a function called "update" and set in a proper maner
    if (node.parent == null) { this.putScatterPoints(node); this.rootnode.add(node); };
    node.visible = true;
}

// loads json data into the scene graph
//      instrdata    -  json instrument definition data
//      particledata  -  json particle ray data
//      main         -  reference to scene wrapper
var TraceLoader = function(instrdata, particledata, main)
{
    this.main = main;
    this.instrdata = instrdata;
    this.particledata = particledata;
}
// synchronously load instrument component draw data
//
TraceLoader.prototype.loadInstr = function()
{
    var main = this.main;

    // utility function
    var drawFunc = function(call, parentnode, color)
    {
        key = call['key'];
        args = call['args'];

        if (key === 'multiline') {
            main.addMultiLine(args, parentnode, color);
        }
        if (key === 'line') {
            main.addMultiLine(args, parentnode, color);
        }
        if (key === 'circle') {
            main.addCircle(args[0], args[1], args[2], args[3], args[4], parentnode, color);
        }
        if (key === 'sphere') {
            main.addSphere(args[0], args[1], args[2], args[3], 32,32, parentnode, color);
        }
        if (key === 'cone') {
            main.addCone(args[0], args[1], args[2], args[3], args[4], args[5], args[6], args[7], 32, parentnode, color);
        }
        if (key === 'cylinder') {
            main.addCylinder(args[0], args[1], args[2], args[3],args[4], args[5], args[6], args[7], args[8], 32, parentnode, color);
        }
        if (key === 'disc') {
            main.addDisc(args[0], args[1], args[2], args[3],args[4], args[5], args[6], 32, parentnode, color);
        }
        if (key === 'annulus') {
            main.addAnnulus(args[0], args[1], args[2], args[3],args[4], args[5], args[6], args[7], 32, parentnode, color);
        }
        if (key === 'new_circle') {
            main.addNewCircle(args[0], args[1], args[2], args[3],args[4], args[5], args[6], 32, parentnode, color);
        }
        if (key === 'box') {
            main.addBox(args[0], args[1], args[2], args[3], args[4], args[5], args[6], parentnode, color);
        }
        if (key === 'polygon') {
            main.addPolygon(args[0], args[1], parentnode, color);
        }
    }

    // INSTRUMENT
    var instr = this.instrdata;
    var instname = instr['name'];
    var abspath = instr['abspath'];

    // COMPONENTS
    var comps = instr['components'];
    var comp;
    var comp_node;
    var compname;
    var m4;
    var acolor;

    for (var i = 0; i < comps.length; i++)
    {
        comp = comps[i];

        compname = comp['name'];
        m4 = comp['m4'];
        comp_node = main.addComponent(compname);
        acolor = main.getNextComponentColor();

        drawcalls = comp['drawcalls']
        var call;
        for (var j = 0; j < drawcalls.length; j++) {
            call = drawcalls[j];

            drawFunc(call, comp_node, acolor);
        }

        comp_matrix = new THREE.Matrix4();
        comp_matrix.set(m4[0], m4[1], m4[2], m4[3], m4[4], m4[5], m4[6], m4[7], m4[8], m4[9], m4[10], m4[11], m4[12], m4[13], m4[14], m4[15]);
        comp_node.applyMatrix4(comp_matrix);
    }
}
// load particle rays
//
TraceLoader.prototype.loadParticles = function()
{
    var main = this.main;

    // set Lut stuff
    main.setLutRange(this.particledata['vmin'], this.particledata['vmax']);

    // RAYS
    var rays = this.particledata['rays'];
    var ray;
    var aVertices;
    for (var i = 0; i < rays.length; i++) {
        ray = rays[i];

        rayobj = new THREE.Object3D();
        aVertices = [];
        var aCompVertices;
        var compname;
        var speed = ray['speed'];

        var groups =  ray['groups'];
        var group;
        for (var j = 0; j < groups.length; j++) {
            group = groups[j];
            compname = group['compname'];
            aCompVertices = [];

            // PARTICLE STATES
            var events = group['events'];
            var localstate;
            var args;
            for (var k = 0; k < events.length; k++) {
                localstate = events[k];

                args = localstate['args'];
                aCompVertices.push(new THREE.Vector3(args[0], args[1], args[2]));
            }

            // transform these vertices by component matrix and add to vertex container for this ray
            if (main.compnodes[compname])
            {
                aVertices = aVertices.concat(transformPoints(aCompVertices, main.compnodes[compname].matrix));
            }
        }
        // add ray as a multiline
        main.addRayNode(rayobj, aVertices, speed);
    }
}
//  program controller
//      campos_x/y/z  -  determines initial camera position, this is used with --inspect
var Controller = function(campos_x, campos_y, campos_z, box_lst, invert_canvas)
{
    this.camPosInitial = new THREE.Vector3(campos_x, campos_y, campos_z);
    this.invert_canvas = invert_canvas;
    this.box_lst = box_lst; // this would be [x_min, x_max, ...]
    this.main = new Main();
    this.loader = new TraceLoader(MCDATA_instrdata, MCDATA_particledata, this.main);
    this.numRays = MCDATA_particledata["numrays"]
    this.viewmodel = new ViewModel(numRays = this.numRays);
}
Controller.prototype.setUpdateGuiFunc = function(updateGuiFunc)
{
    this.updateGuiFunc = updateGuiFunc;
}
//  main program execution loops are set up here
//
Controller.prototype.run = function()
{
    // init mcdisplay
    this.main.init(this.camPosInitial, this.invert_canvas);

    // execution loops
    var _this = this;
    var renderLoop = function()
    {
    	requestAnimationFrame(renderLoop);
    	_this.main.renderer.render(_this.main.scene, _this.main.camera);
    }
    var busy = false;
    var dataLoop = function()
    {
        // add a single-threaded critical section / lock
        if (busy == true) { console.log("busy"); return; }
        busy = true;

        // show/hide ray nodes
        if (_this.viewmodel.playBack == PlayBack.ALL)
        {
            _this.showAllRays();
            return;
        }
        if (_this.viewmodel.playBack == PlayBack.RUN)
        {
            _this.incSingleRay();
        }
        if (_this.viewmodel.playBack == PlayBack.PAUSE)
        {
            // TODO: use viewmodel.updateVersion to skip this step if needed
            _this.showCurrentRay();
        }
        if (viewmodel.getDisplayMode() == DisplayMode.SINGLE)
        {
            _this.hidePrevRays();
        }

        // set scatter cubes on/off
        if (viewmodel.getShowScatterCubes())
        {
            _this.main.showCubes();
        }
        else
        {
            _this.main.hideCubes();
        }

        setTimeout(dataLoop, 1000/_this.viewmodel.raysPrSec);

        busy = false;
    }
    var updateGuiLoop = function()
    {
        setTimeout(updateGuiLoop, 100);
        _this.updateGuiFunc();
    }

    // initiate timed execution loops
    renderLoop();
    dataLoop();
    updateGuiLoop();

    // load data - possibly heavy
    this.loader.loadInstr();
    this.loader.loadParticles();

    // set pause playback mode now after everything else, in case of zero or one rays
    if (this.numRays <= 1) { this.viewmodel.playBack = PlayBack.PAUSE }
}
Controller.prototype.showAllRays = function()
{
    this.main.showAllRays();
}
Controller.prototype.hidePrevRays = function()
{
    var arrLast = this.viewmodel.shiftNonFirstRayIdxs();
    for (var i = 0; i < arrLast.length; i++)
    {
        this.main.hideRay(arrLast[i]);
    }
}
Controller.prototype.incSingleRay = function()
{
    this.viewmodel.setRayIdx(this.viewmodel.getRayIdx() + 1)
    var retidx = this.main.showRay(this.viewmodel.getRayIdx());
}
Controller.prototype.showCurrentRay = function()
{
    this.main.showRay(this.viewmodel.getRayIdx());
}
Controller.prototype.setViewTop = function()
{
    var box = this.box_lst;
    var x = - 0.01; // must not be zero due to the euler angle rotation lock stuff
    var y = (box[5] - box[4])/2;
    var z = (box[5] - box[4])/2;
    this.main.camera.up = new THREE.Vector3(1, 0, 0);
    this.main.setCameraView(new THREE.Vector3(x, y, z));
    this.main.camera.up = new THREE.Vector3(0, 1, 0);
    // next line just to be safe
    this.main.camera.updateMatrix();
}
Controller.prototype.setViewSide = function()
{
    var box = this.box_lst;
    var x = - (box[5] - box[4])/2;
    var y = 0;
    var z = (box[5] - box[4])/2;
    this.main.setCameraView(new THREE.Vector3(x, y, z));
    // next line just to be safe
    this.main.camera.updateMatrix();
}
Controller.prototype.setViewHome = function()
{
    this.main.setCameraView(this.camPosInitial);
}

//  enum for playback state
//
PlayBack = { RUN : 0, PAUSE : 1, ALL : 3 };

//  displaymode e.g. keep rays animated or display one at the time
//
DisplayMode = { SINGLE : 0, KEEP : 1 }

//  viewmodel for keeping control data free of the gui
//
var ViewModel = function(numRays)
{
    this.playBack = PlayBack.RUN;
    this.displayMode = DisplayMode.SINGLE;

    this.showScattCubes = false;

    this.numRays = numRays;
    this.rayIdx = [-1];
    this.raysPrSec = 5;

    this.updateVersion = 0; // incremented on update
}
ViewModel.prototype.getUpdateVersion = function()
{
    return this.updateVersion;
}
ViewModel.prototype.setShowScatterCubes = function(scOnOff)
{
    if (scOnOff)
    {
        this.showScattCubes = true;
    }
    else
    {
        this.showScattCubes = false;
    }
    this.updateVersion += 1;
}
ViewModel.prototype.getShowScatterCubes = function()
{
    return this.showScattCubes;
}
ViewModel.prototype.setRayIdx = function(idx)
{
    // the socalled "js modulus bug" for negative numbers, using ((n % m) + m) % m;
    idx = ((idx % this.numRays) + this.numRays) % this.numRays;
    this.rayIdx.push(idx);
    this.updateVersion += 1;
}
ViewModel.prototype.getRayIdx = function(idx)
{
    return this.rayIdx[this.rayIdx.length-1];
}
ViewModel.prototype.shiftNonFirstRayIdxs = function()
{
    if (this.rayIdx.length <= 1)
    {
        return [];
    }

    var arr = [];
    for (var i=0; i < this.rayIdx.length -1; i++)
    {
        arr.push(this.rayIdx.shift());
    }
    this.updateVersion += 1;
    return arr;
}
ViewModel.prototype.setPlayBack = function(playBack)
{
    this.playBack = playBack;
    this.updateVersion += 1;
}
ViewModel.prototype.getPlayBack = function()
{
    return this.playBack;
}
ViewModel.prototype.setDisplayMode = function(displayMode)
{
    this.displayMode = displayMode;
    this.updateVersion += 1;
}
ViewModel.prototype.getDisplayMode = function()
{
    return this.displayMode;
}

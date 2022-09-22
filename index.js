const express = require('express');
const spawn = require("child_process").spawn;
const app = express();
const port = 3000;
const { v4: uuidv4 } = require('uuid');
const multer = require('multer');
const path = require('path');
const fs = require('fs');
const sqlite3 = require('sqlite3').verbose();
const zmq = require('zeromq');

const db = new sqlite3.Database('./db/SCV.db');

app.use(express.static('static'));
app.use(express.static('public'));

app.listen(port, () => {
    console.log("GLmol listening at http://localhost:"+ port);
});

// const gen3dPython= spawn('python3', ["./python/gen_3d.py"]);
// gen3dPython.stdout.on('data', (data) => {
//     console.log(`gen_3d:\n${data}`);
// });
//
// const genDictPython = spawn('python3', ['./python/gen_dict.py']);
// genDictPython.stdout.on('data', (data) => {
//     console.log(`gen_dict:\n${data}`);
// });


const storage = multer.diskStorage({
    destination: (req, file, cb) => {
        const dir = 'upload/' + req.body.job_number;
        if(!fs.existsSync(dir)){
            return fs.mkdir(dir, error=> cb(error, dir));
        }
        return cb(null, dir);
    },
    filename: (req, file, cb) => {
        console.log(file);
        cb(null, file.originalname);
    }
});
const upload = multer({storage});

const model_requester = zmq.socket('req');
const dict_requester = zmq.socket("req");

model_requester.connect("ipc:///tmp/gen3d.ipc");
dict_requester.connect("ipc:///tmp/genDicts.ipc");

let queue = {}

model_requester.on("message", reply => {
    // console.log(reply);
    if(reply.toString() == "error")
    {
        console.error("error in model requester");
    }
    else {
        queue[JSON.parse(reply)['job_number']].send(reply);
        delete queue[reply['job_number']];
    }
});

app.use(function (err, req, res, next) {
    console.log('This is the invalid field ->', err.field)
    next(err)
});

app.post('/job', upload.array('user_upload', 100), (req, res) => {
    let pdb_dest = "";
    if(req.files.length != 0) {
        pdb_dest = req.files[0].destination;
    }
    let job_number = req.body.job_number;
    let psms = req.body.psms;
    let ptms = req.body.ptms;
    let background_color = req.body.background_color;
    let species = req.body.species;
    console.log(pdb_dest);
    let req_dict = {
        'job_number': job_number,
        'psms': psms,
        'ptms': ptms,
        'background_color': background_color,
        'pdb_dest': pdb_dest,
        'species': species
    }

    dict_requester.send(JSON.stringify(req_dict));
    res.end(job_number);
});

app.get('/view', (req, res) => {
    console.log(req.query);
    res.sendFile(path.join(__dirname, '/static/view.html'));
});

app.post('/protein-list',  upload.none(), (req, res) => {
    let job_number = req.body.job;
    db.get(`SELECT * FROM results WHERE job_number = ?`, [job_number], (err, row) => {
        if (err) {
            console.error(err);
            res.end();
        }
        res.send(JSON.stringify(row));
    })
});


app.post('/3d-view', upload.none(), (req, res) => {
    let job_number = req.body.job;
    let percent = req.body.percent;
    let protein = req.body.protein_id;
    let freqArr = req.body.freqArr;
    let id_ptm_idx_dict = req.body.id_ptm_idx_dict;
    let regex_dict = req.body.regex_dict;
    let background_color = req.body.background_color;
    let pdb_dest = req.body.pdb_dest;
    console.log(background_color);
    console.log(id_ptm_idx_dict);
    let res_path = path.join(__dirname, '/results/' + job_number + '_' + protein + '.json');

    if(fs.existsSync(res_path)){
        res.sendFile(res_path);
    }
    else {
        let request_dict = {
            'job_number': job_number,
            'percent': percent,
            'protein': protein,
            'freqArr': freqArr,
            'id_ptm_idx_dict': id_ptm_idx_dict,
            'regex_dict': regex_dict,
            'background_color': background_color,
            'pdb_dest': pdb_dest
        }
        queue[job_number] = res;
        model_requester.send(JSON.stringify(request_dict));
    }
});
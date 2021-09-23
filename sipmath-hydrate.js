const jStat = require("jstat")
const Library = require('./siplibs/test.json');
// const {GoogleSpreadsheet} = require('google-spreadsheet');
// const fetch = require("node-fetch");
// require('dotenv').config();
// config go to file later

// const for now
const trials = 100;

// start hydrating a sip from the library
async function setupFeeds() {
    console.log("Starting Hydrate")
    //console.log(Library)
    simulateSIP(Library, "thetatoken") // see siplibs/test.json for token names
}

function simulateSIP(selfIn, sip) {
    // Expects library as input and the name of A sip
    // TODO: Add an all option for doing all sips
    // example keyword argument for input: simulateSIP("Variable1", HDR2= {"seed3":0, "seed4":1})

    let randomarray = [];
    let returnValue = [];
    // console.log("sip requested ", sip);
    let sipIndex = selfIn.sips.findIndex((item) => item.name === sip);
    console.log("sipIndex ", sipIndex);

    let aCoeffs = selfIn.sips[sipIndex].arguments.aCoefficients;
    console.log("The aCoeffs for ", aCoeffs);
    let lowerBound = "";
    let upperBound = "";
    let a = "";
    console.log(selfIn.sips[sipIndex].ref);
    let functionName = selfIn.sips[sipIndex].function;

    if (selfIn.sips[sipIndex].ref.source == "copula") {
        console.log("matching copula");
        randomarray = generateCopula(selfIn, selfIn.sips[sipIndex].ref.copulaLayer); // c1 or c2 etc
        console.log("copula generate ", randomarray);
    } else if (selfIn.sips[sipIndex].ref.source == "rng") {
        console.log("hmm matching rng sip? that's ok");
        // need "hdr1" and library as input to generateRandom
         randomarray = generateRandom(
            selfIn.sips[sipIndex].ref.name,
            selfIn.U01
        );
        //console.log("yo dog ", randomarray);
        // randomarray[0] = temprandomarray;
        console.log("NOT copula sip ", randomarray);
    }

    try {
        // Something that could throw exception TODO avoid this
        lowerBound = selfIn.sips[sipIndex]["arguments"]["lowerBound"];
    } catch (e) {
        console.log("Nothing to see here. Just no lowerBound");
    }
    try {
        upperBound = selfIn.sips[sipIndex]["arguments"]["upperBound"];
    } catch (e) {
        console.log("Nothing to see here. Just no lowerBound");
    }

    //let functionName = selfIn.sips[sipIndex].function;

    // if the function is a built in function Metalog_1_0 then,
    // for each item
    if (functionName == "Metalog_1_0") {
        // TODO: change for loop into a vectorized numpy function
        if (randomarray.constructor === Array) { //hmm forgot why I did this. I'm getting old
            let wrapArrayInArray = [randomarray]; // looking for array of arrays
            console.log("wrapArrayInArray ", wrapArrayInArray);
            console.log("process non copula sip");   
            // console.log("randomarray ", randomarray);
            for (var i = 0; i < randomarray.length; i++) {
                console.log("randomarray[i] ", randomarray[i]);
                let ml = metalog(randomarray[i], aCoeffs, lowerBound, upperBound);                //console.log("a ", a);
                returnValue.push(ml);
            }

        }
        console.log("sipSim length of randomarray ", randomarray.length);
        for (let index = 0; index < randomarray.length; ++index) {
            let ml = metalog(randomarray[index], aCoeffs, lowerBound, upperBound);
            //returnValue.append(ml);
            returnValue.push(ml);
        }
        //for trial in randomarray:
    }
    //selfIn.sips = returnValue;
    //return randomarray;
    return returnValue;
}


function metalog(y, a, bl = "", bu = "") {
    // y = array of unis, a = aacoeffs
    //let t = Array(a.length); // a and now t are acoeff
    console.log("how many of aCoeffs ", a.length);

    function convert_to_float(a) {
        // Using parseFloat() method
        var floatValue = parseFloat(a);
        // Return float value
        return floatValue;
    }

    //Array(6).fill(null).map((_, i) => i);
    console.log("y ", typeof y);
    console.log("y me ", y); // y[x] is an array of unis may have been thru copula calc
    let vector = [];
    // np_a = np.array(a).reshape(-1, 1)
    console.log("a ", typeof a, "a length ", a.length);

    console.log("a ", a);
    let np_a = a; // !!TODO CONFIRM matrixmult does this need to be vertical array??
    //for n in t:
    for (let index = 1; index < (a.length + 1); ++index) { //cant start with 0 coeeffs for each aCoeff
                vector.push(basis(y, index));
        }    

    //vector = np.array(vector, (dtype = object));

    // let mky = np.matmul(vector, np_a);
    console.log("np_a ", typeof np_a, 'length: ', np_a.length);
    console.table( np_a);
    console.log("vector ", typeof vector, " length ", vector.length);
    console.table( vector);

/*     let MatrixProd = (A, B) =>
    A.map((row, i) =>
      B[0].map((_, j) =>
        row.reduce((acc, _, n) =>
          acc + A[i][n] * B[n][j], 0
        )
      )
    ) */
//let A = [[8, 3], [2, 4], [3, 6]];
//let B = [[1, 2, 3], [4, 6, 8]];
console.log('this is B ', typeof np_a, np_a);
//let mky = MatrixProd(vector, np_a);
// np_a = [ [1], [-3.37926448792588], [1.5782860765295803], [-0.4670501767969924] ]
//vector =[[0.994838102415617],[0.004685644724859045],[-0.004451370153069208],[-0.013144268586090958]]
//above seems to work
    let wrappedVector = [vector];
    let wrappedNp_x = [np_a];
    let herItis = [[]];
    let wrappedNp_a = wrappedNp_x[0].map(e => [e])

    /* wrappedNp_a.forEach(occur => {
        wrappedNp_a[occur] = occur.map(convert_to_float);
    }); */
    //wrappedNp_a.forEach(function(is) { is.forEach( function (his, i) { console.log(his, i) } ) } );
    console.log("wrappedVector ", wrappedVector);
    console.log("wrappedNp_a ", wrappedNp_a);
    let mky = multiply(wrappedVector, wrappedNp_a);
    console.log("mm array of unis with acoeffs ", mky);

    // Unbounded
    // if (type(bl) == str) and (type(bu) == str):
    if (
        (typeof bl === "string" || bl instanceof String) &&
        (typeof bu === "string" || bu instanceof String)
    ) {
        return mky;
    }
    // Bounded lower
    else if (typeof bl !== "string" && typeof bu == "string") {
        convert_to_float(bl);
        return bl + Math.exp(mky);
    }
    // Bounded upper
    else if (typeof bl == "string" && typeof bu != "string") {
        convert_to_float(bu);
        return bu - Math.exp(-mky);
    }
    // Bounded
    else if (typeof bl != "string" && typeof bu != "string") {
        return bl + (bu * Math.exp(mky)) / (1 + Math.exp(mky));
    }
}
function generateCopula(selfy, copulaCount) {
    console.log(selfy, copulaCount);
    let ret = [];
    let whichCorrelationMatrix = [];
    // Do all of this for all copulas in document
    selfy.U01.copula.forEach((copula) => {
        if (copula.function == "GaussianCopula") {
            // now get the cholesky factors
            //  get the global variable
            let specifyCorrelationMatrix = copula.arguments.correlationMatrix.value;
            let copulaArgs = copula.arguments.rng;
            let randomMatrixOfHDRs = [];
            for (let i = 0; i < copulaArgs.length; i++) {
                console.log("caller ", copulaArgs[i]);
                let val = generateRandom(copulaArgs[i], selfy.U01); // from U01/RNG
                /*  {
                            "counter": "PM_Index",
                            "entity": 1,
                            "varId": 6187319,
                            "seed3": 0,
                            "seed4": 0
                        } */

                randomMatrixOfHDRs.push(val);
                console.log("randomMatrixOfHDRs ", randomMatrixOfHDRs);
            }
            // WORKS as of here
            selfy.globalVariables.forEach((item, index) => {
                if (item["name"] == specifyCorrelationMatrix) {
                    whichCorrelationMatrix = index;
                } else {
                    let index = -1;
                }
            });

            let thisCorrelationMatrix =
                selfy.globalVariables[whichCorrelationMatrix].value.matrix;
            let correlationMatrix = convertMx(thisCorrelationMatrix);

            // Find the Cholesky Factors
            //let lu = jStat.lu(correlationMatrix);
            // console.log("you are lower tri ", correlationMatrix);
            let cho = jStat(jStat.cholesky(correlationMatrix));

            //cho = cholesky(correlationMatrix, (lower = False)); //  TODO: jstat?
            console.log("Cholesky: \n", cho);
            // LOOKS GOOD HERE BUT NEED TO CONFIRM

            //cho = np.matrix(cho)
            //Apply the Cholesky Factors to the randoms
            let col = copula.copulaLayer.findIndex((item) => item === copulaCount); //
            console.log("this is copula level eg 2=c3,3=c4 ", col);
            let choSubSample = cho[col].slice(0, col + 1); //
            console.log("this is cho sample ", choSubSample);
            let runiRow = [];
            let corrSamples = [];

            for (let i = 0; i < randomMatrixOfHDRs[0].length; i++) {
                // for each trail
                let randomMatrixHRDsSample = [];
                for (let index = 0; index < col + 1; ++index) {
                    //each variable upto pos col
                    // get first x cols in randuniframe
                    randomMatrixHRDsSample[index] = randomMatrixOfHDRs[index]; //
                    //console.log(randomMatrixHRDsSample[index].length);
                    // console.log("runi samples only arrays needed ", randomMatrixHRDsSample);

                    runiRow[i] = randomMatrixHRDsSample.map(function (x) {
                        //console.log("runi row ", index, x[0]);
                        return x[i];
                    });
                    //          console.log("runi row ", runiRow);
                }
                // console.log("runi row and choSubSample ", runiRow[i], choSubSample);

                // GOOD TO HERE?
                corrSamples[i] = jStat.dot(runiRow[i], choSubSample);
                // let mMult = jStat.dot(invCdf, choSubSample);
            }
            console.log("YIPPE check me... ", corrSamples);
            ret = corrSamples;
            //for (let i = 0; i < randomMatrixHRDsSample[index].length; i++) {
            //  corrSamples[index] = cho.multiply(
            //    randomMatrixHRDsSample[index][i],
            //    choSubSample
            //  );
            //console.log("corrSamples ", corrSamples);
            // changed to trials ?? why 1000?? hmm
            // asked for
            // if C4 we need uniRandoArrays[] 0,1,2,3
            // find out what rows of the hdr ret matrix to use based on what was asked for
            //           col = copulas["copulaLayer"].index(random)
            //           choSubSample = cho[:col+1, col]
            //           trial = randomMatrix[1:col+2, i]
            //           invCdf = norm.ppf(trial).reshape(-1, col+1)
            //           # get HDRs from matrix
            //           mMult = np.dot(invCdf, choSubSample)
            //           val = float(norm.cdf(mMult))
            //           ret.append(val)

            //console.log(randomMatrixOfHDRs);
            //console.log("matrix mult ", corrSamples);
            //let choSubSample = cho[(0, [col + 1, col])];
            //console.log("this is chosample ", choSubSample);
            //  let trial = randomMatrixOfHDRs[(1, [col + 2, i])];

            // let invCdf = norm.ppf(trial).reshape(-1, col + 1);
            // get HDRs from matrix
            // let mMult = jStat.dot(invCdf, choSubSample);
            // let val = norm.cdf(mMult);

            //Convert to uniform correlated samples over 0,1 using normal CDF
            // var uniformCorrSamples = corrSamples.map(function (x) {
            // var normDist = jStat.normal(0, 1); //FTW Cholesky needs normal?
            //   var normDist = jStat.normal.inv(x, 0, 1);
            //  hrd , 0,1 --- If not cormatrix - but rnd straight into metalog else output of cho goes into metalog
            // return normDist.cdf(x);
            //console.log(normDist);
            //   return normDist;
            //   });
            // }
            // }
            //  let val = uniformCorrSamples;
            //console.log(val);
            //  ret.push(val);
            //  }
        } else {
            console.log(
                "TypeError The function type for this copula is unsupported. "
            );
        }
    });
    // Add the result of the cholesky to the copula variable ??
    //self.copula = corrSamples;
    return ret;
}

function generateRandom(args, selfIn) {
    // from U01/RNG
    /*  {
                      "counter": "PM_Index",
                      "entity": 1,
                      "varId": 6187319,
                      "seed3": 0,
                      "seed4": 0
                  } */
    //console.log("in bound ",args);
    // console.log("in bound selfIn ", selfIn.rng);
    let rngArgs = selfIn.rng.findIndex((x) => x.name === args);
    // console.log(args, rngArgs);
    console.log("seed ", selfIn.rng[rngArgs].arguments.varId);

    //let rngArgs = self.U01.rng.indexOf(name.args);
    //console.log("use as seed ", rngArgs.arguments); // use as seeds
    var samples = [];
    const seedPerDist = selfIn.rng[rngArgs].arguments.varId;
    //var seedPerDist = Math.random();
    // console.log(seedPerDist); // just one for now
    for (var distTrials = 0; distTrials < trials; distTrials++) {
        samples[distTrials] = HDRando(seedPerDist, distTrials);
    }
    console.log("samples og ", samples);
    return samples;
}

// hubbardresearch.com for more info
function HDRando(seed, PM_Index) {
    const largePrime = 2147483647;
    const million = 1000000;
    const tenMillion = 10000000;
    // Do we need this in js? or is there a modulo?
    function MOD(n, m) {
        var remain = n % m;
        return Math.floor(remain >= 0 ? remain : remain + m);
    }
    let randi =
        (MOD(
                (MOD(
                        (seed + million) ^ (2 + (seed + million) * (PM_Index + tenMillion)),
                        99999989
                    ) +
                    1000007) *
                (MOD(
                        (PM_Index + tenMillion) ^
                        (2 +
                            (PM_Index + tenMillion) *
                            MOD(
                                (seed + million) ^
                                (2 + (seed + million) * (PM_Index + tenMillion)),
                                99999989
                            )),
                        99999989
                    ) +
                    1000013),
                largePrime
            ) +
            0.5) /
        largePrime;
    return randi;
}

function convertMx(correlationMatrix) {
    let variables = [];

    //gotta figure out all of the variables in the matrix
    //for vars in correlationMatrix:
    //    if vars["row"] not in variables:
    //        variables.append(vars["row"])
    //console.log("corm ", correlationMatrix);
    correlationMatrix.forEach((sipVar) => {
        //console.log("here too ", sipVar);
        if (variables.includes(sipVar["row"])) {
            //console.log("variables already in");
        } else {
            variables.push(sipVar["row"]);
            //console.log("here now ", variables);
        }
    });

    let variableCount = variables.length;

    // let returnArray = np.zeros(shape=[variableCount, variableCount])
    let returnArray = Array(variableCount)
        .fill()
        .map(() => Array(variableCount).fill(0));
    //console.log("returnarr", returnArray);

    correlationMatrix.forEach((items) => {
        //console.log("vars ", variables);
        //console.log("items ", items);

        let i = items.row;
        let j = items.col;

        i = variables.indexOf(items["row"]);
        j = variables.indexOf(items["col"]);
        let value = items.value;
        returnArray[i][j] = value;
        returnArray[j][i] = value;
    });
    /*
    var testMatrix = [
      [1, 0.847718, 0.055669605, 0.485944056, 0.564808183],
      [0, 0.530447163, 0.029540807, 0.247911917, 0.045299211],
      [0, 0, 0.998012142, -0.28349892, -0.080721634],
      [0, 0, 0, 0.788686514, 0.186186642],
      [0, 0, 0, 0, 0.798597678]
    ];
    */
    //console.log("ret array ", returnArray);
    return returnArray;
}

function multiply(a, b) {
    console.log("a ", a, "b",b);
    var aNumRows = a.length,
        aNumCols = a[0].length || 0, // if a is a vector
        bNumRows = b.length,
        bNumCols = b[0].length || 0,
        m = new Array(aNumRows); // initialize array of rows
    console.log("aNumRows ", aNumRows, "aNumCols ", aNumCols,"bNumCols",bNumCols,"bNumCols",bNumCols);
    for (var r = 0; r < aNumRows; ++r) {
        m[r] = new Array(bNumCols); // initialize the current row
        for (var c = 0; c < bNumCols; ++c) {
            console.log("r ", r, "c ", c);
            m[r][c] = 0; // initialize the current cell
            for (var i = 0; i < aNumCols; ++i) {
                m[r][c] += a[r][i] * b[i][c];
                console.log("a[r][i] ", a[r][i]);
            }
        }
    }
    return m;
}

function basis(y, t) {
    console.log("aCoeff position in basis ", t, y);
    //console.log("y in basis ", y);
    let ret = 0;
    if (t == 1) {
        ret = 1;
    } else if (t == 2) {
        ret = Math.log(y / (1 - y));
        // console.log("ret when t2 ", y, ret);
        if (isNaN(ret)) {
            console.log("ret when t2 ", y, ret);
        }
    } else if (t == 3) {
        ret = (y - 0.5) * Math.log(y / (1 - y));
        if (isNaN(ret)) {
            console.log("ret when t3 ", y, ret);
        }
    } else if (t == 4) {
        ret = y - 0.5;
        if (isNaN(ret)) {
            console.log("ret when t4 ", y, ret);
        }
    } else if (t >= 5 && t % 2 == 1) {
        ret = Math.pow(y - 0.5, Math.floor((t - 1) / 2));
        if (isNaN(ret)) {
            console.log("ret when t5 ", y, ret);
        }
    } // requires js int division
    else if (t >= 6 && t % 2 == 0) {
        ret = Math.pow(y - 0.5, Math.floor((t - 1) / 2)) * Math.log(y / (1 - y));
        if (isNaN(ret)) {
            console.log("ret when t6>= ", y, ret);
        }
    } // requires js int division
    console.log("out from basis ", ret);
    return ret;
}

//setupContract()
setupFeeds()
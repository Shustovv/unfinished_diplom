  window.addEventListener('load', () => {
    const canvasElements = document.querySelectorAll('#canvas-container canvas');
    const canvas = document.getElementById('canvas1');
    const ctx = canvas.getContext('2d');
    const canvas2 = document.getElementById('canvas2');
    const ctx2 = canvas2.getContext('2d');
    const canvas3 = document.getElementById('canvas3');
    const ctx3 = canvas3.getContext('2d');

    function showCanvas(canvasNumber) {
      // Скрыть все canvas
      canvasElements.forEach(canvas => {
        canvas.style.display = 'none';
      });
      // Показать выбранную canvas
      const selectedCanvas = document.getElementById(`canvas${canvasNumber}`);
      selectedCanvas.style.display = 'block';
    }

    // Обработчик событий для ссылок в sidebar
    const sidebarLinks = document.querySelectorAll('.canvas-list a');
    sidebarLinks.forEach((link, index) => {
      link.addEventListener('click', (event) => {
        event.preventDefault();
        showCanvas(index + 1);
      });
    });
      
    function buildBeam() {
      
      const width = parseFloat(document.getElementById('width').value);
      const height = parseFloat(document.getElementById('height').value);
      const xDivisions = parseInt(document.getElementById('x-divisions').value);
      const yDivisions = parseInt(document.getElementById('y-divisions').value);
      const loadAreaX1 = parseInt(document.getElementById('load-area-x1').value);
      const loadAreaX2 = parseInt(document.getElementById('load-area-x2').value);
      const beamFixing = document.getElementById('beam-fixing').value;
      const beamType = document.getElementById('beam-type').value;
      const poisson = parseFloat(document.getElementById('material-poisson-ratio').value);
      const ung = parseFloat(document.getElementById('material-modulus').value);
      const force = parseInt(document.getElementById('load-force').value);
      const xStep = width / xDivisions;
      const yStep = height / yDivisions;
      const nodes = [];
      for (let i = 0; i <= yDivisions; i++) {
        for (let j = 0; j <= xDivisions; j++) {
          const x = j * xStep;
          const y = i * yStep;
          nodes.push({ x, y });
        }
      } 
      
      const stress = [
         0.0133308,
         0.0156909,
         0.0366498,
         0.0690677,
         0.120336,
         0.197416,
         0.301598,
         0.42049,
         0.534735,
         0.623613,
         0.670638,
         0.670611,
         0.623517,
         0.534536,
         0.420103,
         0.300869,
         0.196051,
         0.117794,
         0.0643664,
         0.0280119,
         0.024103,
         0.0261087,
         0.0307666,
         0.0719717,
         0.135984,
         0.237788,
         0.391888,
         0.601602,
         0.841531,
         1.07232,
         1.25189,
         1.34659,
         1.34654,
         1.2517,
         1.07193,
         0.840776,
         0.600178,
         0.389224,
         0.232828,
         0.126796,
         0.0550669,
         0.0472989,
         0.0500264,
         0.0590942,
         0.138658,
         0.26334,
         0.463896,
         0.771707,
         1.19665,
         1.68523,
         2.15615,
         2.52267,
         2.71461,
         2.71451,
         2.52231,
         2.15541,
         1.68379,
         1.19393,
         0.766631,
         0.454437,
         0.245793,
         0.106301,
         0.0909875,
         0.0696686,
         0.082629,
         0.194868,
         0.373339,
         0.666025,
         1.12616,
         1.77785,
         2.53312,
         3.26346,
         3.83223,
         4.12611,
         4.12596,
         3.83175,
         3.26243,
         2.53113,
         1.77411,
         1.11916,
         0.652953,
         0.34903,
         0.149877,
         0.127539,
         0.0831809,
         0.0992036,
         0.235634,
         0.457105,
         0.830634,
         1.4393,
         2.33645,
         3.38693,
         4.40711,
         5.20266,
         5.60431,
         5.60414,
         5.20209,
         4.40591,
         3.3846,
         2.33205,
         1.43104,
         0.81519,
         0.42829,
         0.182033,
         0.153633,
         0.0890794,
         0.106977,
         0.256449,
         0.505732,
         0.942364,
         1.69058,
         2.86136,
         4.24798,
         5.6014,
         6.65978,
         7.17405,
         7.17387,
         6.65918,
         5.60013,
         4.24552,
         2.85672,
         1.68187,
         0.926029,
         0.475138,
         0.199195,
         0.166339,
         0.0864213,
         0.104625,
         0.25364,
         0.510652,
         0.983554,
         1.85205,
         3.33841,
         5.11697,
         6.8625,
         8.23556,
         8.86176,
         8.86159,
         8.23499,
         6.8613,
         5.11463,
         3.33399,
         1.84373,
         0.96792,
         0.481251,
         0.198249,
         0.163414,
         0.0749804,
         0.0915845,
         0.224935,
         0.464463,
         0.934128,
         1.88373,
         3.74962,
         5.99393,
         8.20878,
         9.97208,
         10.6948,
         10.6946,
         9.97159,
         8.20775,
         5.99194,
         3.74585,
         1.87663,
         0.920748,
         0.439197,
         0.17701,
         0.143716,
         0.0554406,
         0.0683374,
         0.17028,
         0.362474,
         0.772358,
         1.72594,
         4.07329,
         6.87804,
         9.66039,
         11.9284,
         12.699,
         12.6989,
         11.928,
         9.65965,
         6.87659,
         4.07055,
         1.72077,
         0.762588,
         0.343957,
         0.134932,
         0.107719,
         0.0295297,
         0.0367046,
         0.0928797,
         0.205394,
         0.477742,
         1.28556,
         4.28585,
         7.76763,
         11.2372,
         14.1933,
         14.8935,
         14.8934,
         14.1931,
         11.2368,
         7.76687,
         4.2844,
         1.28283,
         0.472587,
         0.1956,
         0.0740887,
         0.0580789,
         0.0152647,
         0.0191116,
         0.0491341,
         0.113611,
         0.294871,
         0.980439,
         4.36139,
         8.21353,
         12.0592,
         15.4163,
         16.0399,
         16.0399,
         15.4162,
         12.059,
         8.21314,
         4.36065,
         0.979039,
         0.292224,
         0.108573,
         0.0394338,
         0.0303333,
        ];
        
        function getColorFromGradient(value) {
          const colors = [
            [0, 0, 255],  // Самый светлый цвет (напряжение = minStress)
            [255, 0, 0]         // Самый темный цвет (напряжение = maxStress)
          ];
      
          const r = colors[0][0] * (1 - value) + colors[1][0] * value;
          const g = colors[0][1] * (1 - value) + colors[1][1] * value;
          const b = colors[0][2] * (1 - value) + colors[1][2] * value;
      
          return `rgb(${r}, ${g}, ${b})`;
        }

////////////////////////////////////////////////////////////////////////////////////////

const minStress = Math.min(...stress);
const maxStress = Math.max(...stress);

const data = {
  datasets: [
    {
      type: 'scatter',
      data: nodes,
      pointRadius: 10,
      pointBackgroundColor: stress.map(value => {
        const normalizedValue = (value - minStress) / (maxStress - minStress);
        const color = getColorFromGradient(normalizedValue); // Функция, которая возвращает цвет из градиента
        return color;
      })
    }
  ]
};

      for (let i = 0; i < yDivisions+1; i++){ 
        const startIndex = i * (xDivisions + 1);
        const endIndex = (i + 1) * (xDivisions + 1);
        const lineData = nodes.slice(startIndex, endIndex);
        data.datasets.push({
          type: 'line',
          data: lineData,
          borderColor: 'blue',
          borderWidth: 1,
          fill: false
        });
      }
      for (let i = 0; i <= xDivisions; i++) {
        const lineData = [];
        for (let j = i; j < nodes.length; j += (xDivisions + 1)) {
          lineData.push(nodes[j]);
        }
        data.datasets.push({
          type: 'line',
          data: lineData,
          borderColor: 'blue',
          borderWidth: 1,
          fill: false
        });
        //Линия действия силы
        data.datasets.push({
          type: 'line',
          data: [
            { x: loadAreaX1, y: height },
            { x: loadAreaX2, y: height }
          ],
          borderColor: 'red',
          borderWidth: 4,
          fill: false,
          label: false
        });
      }
      const options = {
        scales: {
          x: {
            type: 'linear',
            position: 'bottom',
            ticks: {
              beginAtZero: true
            }
          },
          y: {
            type: 'linear',
            position: 'left',
            ticks: {
              beginAtZero: true
            }
          }
        },
        plugins: {
          legend: {
            display: false // Скрыть легенду
          }
        }
      };
      if (beamFixing === 'both') {
        // Добавляем левую и правую границы
        data.datasets.push({
          type: 'line',
          data: [
            { x: 0, y: 0 },
            { x: 0, y: height }
          ],
          borderColor: 'orange',
          borderWidth: 7,
          fill: false,
          label: 'beam-boundary'
        });
        data.datasets.push({
          type: 'line',
          data: [
            { x: width, y: 0 },
            { x: width, y: height }
          ],
          borderColor: 'orange',
          borderWidth: 7,
          fill: false,
          label: 'beam-boundary'
        });
      } else if (beamFixing === 'right') {
        // Добавляем правую границу
        data.datasets.push({
          type: 'line',
          data: [
            { x: width, y: 0 },
            { x: width, y: height }
          ],
          borderColor: 'orange',
          borderWidth: 7,
          fill: false,
          label: 'beam-boundary'
        });
      } else if (beamFixing === 'left') {
        // Добавляем левую границу
        data.datasets.push({
          type: 'line',
          data: [
            { x: 0, y: 0 },
            { x: 0, y: height }
          ],
          borderColor: 'orange',
          borderWidth: 7,
          fill: false,
          label: 'beam-boundary'
        });
      }
      new Chart(ctx, {
        type: 'line',
        data: data,
        options: options
      });

      
      const points = [
        { x: 0, y:0},
        { x: 0, y:0.332419},
        { x: 0, y:0.664134},
        { x: 0, y:0.994184},
        { x: 0, y:1.32104},
        { x: 0, y:1.64219},
        { x: 0, y:1.95358},
        { x: 0, y:2.24927},
        { x: 0, y:2.52163},
        { x: 0, y:2.76328},
        { x: 0, y:2.96897},
        { x: 0, y:3.13607},
        { x: 0, y:3.26459},
        { x: 0, y:3.35713},
        { x: 0, y:3.41895},
        { x: 0, y:3.45742},
        { x: 0, y:3.48012},
        { x: 0, y:3.49295},
        { x: 0, y:3.49988},
        { x: 0, y:3.50324},
        { x: 0, y:3.50423},
        { x: 0, y:0},
        { x: 0, y:0.332504},
        { x: 0, y:0.664332},
        { x: 0, y:0.994563},
        { x: 0, y:1.32172},
        { x: 0, y:1.64333},
        { x: 0, y:1.9554},
        { x: 0, y:2.25194},
        { x: 0, y:2.52512},
        { x: 0, y:2.76737},
        { x: 0, y:2.97335},
        { x: 0, y:3.14046},
        { x: 0, y:3.26868},
        { x: 0, y:3.36063},
        { x: 0, y:3.42162},
        { x: 0, y:3.45924},
        { x: 0, y:3.48127},
        { x: 0, y:3.49365},
        { x: 0, y:3.5003},
        { x: 0, y:3.50352},
        { x: 0, y:3.50447},
        { x: 0, y:0},
        { x: 0, y:0.332749},
        { x: 0, y:0.664911},
        { x: 0, y:0.995673},
        { x: 0, y:1.3237},
        { x: 0, y:1.64671},
        { x: 0, y:1.96083},
        { x: 0, y:2.25995},
        { x: 0, y:2.53565},
        { x: 0, y:2.7797},
        { x: 0, y:2.98658},
        { x: 0, y:3.15369},
        { x: 0, y:3.28101},
        { x: 0, y:3.37116},
        { x: 0, y:3.42965},
        { x: 0, y:3.46469},
        { x: 0, y:3.48469},
        { x: 0, y:3.49571},
        { x: 0, y:3.50154},
        { x: 0, y:3.50434},
        { x: 0, y:3.50517},
        { x: 0, y:0},
        { x: 0, y:0.333136},
        { x: 0, y:0.665825},
        { x: 0, y:0.997435},
        { x: 0, y:1.32688},
        { x: 0, y:1.65218},
        { x: 0, y:1.96974},
        { x: 0, y:2.27332},
        { x: 0, y:2.55337},
        { x: 0, y:2.80047},
        { x: 0, y:3.00886},
        { x: 0, y:3.17596},
        { x: 0, y:3.30178},
        { x: 0, y:3.38889},
        { x: 0, y:3.44304},
        { x: 0, y:3.47363},
        { x: 0, y:3.49022},
        { x: 0, y:3.49899},
        { x: 0, y:3.50351},
        { x: 0, y:3.50563},
        { x: 0, y:3.50626},
        { x: 0, y:0},
        { x: 0, y:0.333629},
        { x: 0, y:0.666999},
        { x: 0, y:0.999717},
        { x: 0, y:1.33105},
        { x: 0, y:1.6595},
        { x: 0, y:1.98194},
        { x: 0, y:2.29209},
        { x: 0, y:2.57855},
        { x: 0, y:2.83002},
        { x: 0, y:3.04052},
        { x: 0, y:3.20763},
        { x: 0, y:3.33134},
        { x: 0, y:3.41408},
        { x: 0, y:3.46182},
        { x: 0, y:3.48587},
        { x: 0, y:3.4976},
        { x: 0, y:3.5033},
        { x: 0, y:3.50605},
        { x: 0, y:3.50729},
        { x: 0, y:3.50765},
        { x: 0, y:0},
        { x: 0, y:0.334187},
        { x: 0, y:0.668332},
        { x: 0, y:1.00234},
        { x: 0, y:1.33592},
        { x: 0, y:1.66827},
        { x: 0, y:1.99707},
        { x: 0, y:2.31628},
        { x: 0, y:2.61163},
        { x: 0, y:2.86888},
        { x: 0, y:3.08207},
        { x: 0, y:3.24918},
        { x: 0, y:3.3702},
        { x: 0, y:3.44717},
        { x: 0, y:3.48604},
        { x: 0, y:3.50104},
        { x: 0, y:3.50646},
        { x: 0, y:3.50833},
        { x: 0, y:3.50896},
        { x: 0, y:3.50916},
        { x: 0, y:3.50921},
        { x: 0, y:0},
        { x: 0, y:0.334755},
        { x: 0, y:0.669702},
        { x: 0, y:1.00507},
        { x: 0, y:1.34112},
        { x: 0, y:1.67795},
        { x: 0, y:2.01458},
        { x: 0, y:2.34593},
        { x: 0, y:2.65325},
        { x: 0, y:2.91772},
        { x: 0, y:3.13411},
        { x: 0, y:3.30122},
        { x: 0, y:3.41905},
        { x: 0, y:3.4888},
        { x: 0, y:3.51571},
        { x: 0, y:3.5186},
        { x: 0, y:3.51623},
        { x: 0, y:3.51369},
        { x: 0, y:3.51198},
        { x: 0, y:3.51107},
        { x: 0, y:3.5108},
        { x: 0, y:0},
        { x: 0, y:0.335279},
        { x: 0, y:0.670976},
        { x: 0, y:1.00766},
        { x: 0, y:1.34618},
        { x: 0, y:1.6878},
        { x: 0, y:2.03361},
        { x: 0, y:2.3811},
        { x: 0, y:2.70439},
        { x: 0, y:2.97747},
        { x: 0, y:3.19742},
        { x: 0, y:3.36453},
        { x: 0, y:3.4788},
        { x: 0, y:3.53995},
        { x: 0, y:3.5509},
        { x: 0, y:3.53767},
        { x: 0, y:3.52615},
        { x: 0, y:3.51888},
        { x: 0, y:3.51483},
        { x: 0, y:3.51284},
        { x: 0, y:3.51225},
        { x: 0, y:0},
        { x: 0, y:0.335702},
        { x: 0, y:0.672016},
        { x: 0, y:1.00981},
        { x: 0, y:1.35055},
        { x: 0, y:1.69682},
        { x: 0, y:2.05276},
        { x: 0, y:2.42181},
        { x: 0, y:2.76653},
        { x: 0, y:3.04927},
        { x: 0, y:3.27282},
        { x: 0, y:3.43993},
        { x: 0, y:3.55061},
        { x: 0, y:3.6021},
        { x: 0, y:3.59163},
        { x: 0, y:3.55685},
        { x: 0, y:3.53523},
        { x: 0, y:3.52336},
        { x: 0, y:3.51719},
        { x: 0, y:3.51428},
        { x: 0, y:3.51343},
        { x: 0, y:0},
        { x: 0, y:0.335979},
        { x: 0, y:0.672702},
        { x: 0, y:1.01126},
        { x: 0, y:1.35362},
        { x: 0, y:1.70372},
        { x: 0, y:2.06979},
        { x: 0, y:2.46811},
        { x: 0, y:2.84201},
        { x: 0, y:3.13456},
        { x: 0, y:3.36121},
        { x: 0, y:3.52833},
        { x: 0, y:3.6359},
        { x: 0, y:3.67759},
        { x: 0, y:3.63794},
        { x: 0, y:3.5739},
        { x: 0, y:3.54216},
        { x: 0, y:3.5265},
        { x: 0, y:3.51877},
        { x: 0, y:3.51522},
        { x: 0, y:3.51419},
        { x: 0, y:0},
        { x: 0, y:0.336075},
        { x: 0, y:0.672944},
        { x: 0, y:1.01178},
        { x: 0, y:1.35481},
        { x: 0, y:1.70692},
        { x: 0, y:2.08098},
        { x: 0, y:2.52001},
        { x: 0, y:2.93459},
        { x: 0, y:3.23502},
        { x: 0, y:3.46341},
        { x: 0, y:3.63052},
        { x: 0, y:3.73636},
        { x: 0, y:3.77017},
        { x: 0, y:3.68985},
        { x: 0, y:3.5851},
        { x: 0, y:3.54538},
        { x: 0, y:3.52772},
        { x: 0, y:3.51935},
        { x: 0, y:3.51555},
        { x: 0, y:3.51446},
        ];
   
      //.........................График деформаций..............................................................................
        
      const options2 = {
        scales: {
          x: {
            type: 'linear',
            position: 'bottom',
            ticks: {
              beginAtZero: true
            }
          },
          y: {
            type: 'linear',
            position: 'left',
            ticks: {
              beginAtZero: true
            }
          }
        },
        plugins: {
          legend: {
            display: false // Скрыть легенду
          }
        }
      };
      
        const summedVectors = [];
      
        for (let i = 0; i < nodes.length; i++) {
          const summedVector = {
            x: nodes[i].x - points[i].x,
            y: nodes[i].y - points[i].y
        };
        summedVectors.push(summedVector);
        }
        //строится график деформаций
        const data2 = {
          datasets:[
            {
              type: 'scatter',
              data: summedVectors,
              pointRadius: 3,
              pointBackgroundColor: 'red',
              pointHoverRadius: 8
            }
          ]
        };  
        if (beamFixing === 'both') {
          // Добавляем левую и правую границы
          data2.datasets.push({
            type: 'line',
            data: [
              { x: 0, y: 0 },
              { x: 0, y: height }
            ],
            borderColor: 'orange',
            borderWidth: 7,
            fill: false,
          });
          data2.datasets.push({
            type: 'line',
            data: [
              { x: width, y: 0 },
              { x: width, y: height }
            ],
            borderColor: 'orange',
            borderWidth: 7,
            fill: false,
          });
        } else if (beamFixing === 'right') {
          // Добавляем правую границу
          data2.datasets.push({
            type: 'line',
            data: [
              { x: width, y: 0 },
              { x: width, y: height }
            ],
            borderColor: 'orange',
            borderWidth: 7,
            fill: false,
          });
        } else if (beamFixing === 'left') {
          // Добавляем левую границу
          data2.datasets.push({
            type: 'line',
            data: [
              { x: 0, y: 0 },
              { x: 0, y: height }
            ],
            borderColor: 'orange',
            borderWidth: 7,
            fill: false,
          });
        }
        // Добавляем горизонтальные линии
      for (let i = 0; i < yDivisions+1; i++) {
        startIndex = i * (xDivisions + 1);
        endIndex = (i + 1) * (xDivisions + 1);
        const lineData2 = summedVectors.slice(startIndex, endIndex);
        data2.datasets.push({
          type: 'line',
          data: lineData2,
          borderColor: 'blue',
          borderWidth: 1,
          fill: false
        });
      }
      // Добавляем вертикальные линии
      for (let i = 0; i <= xDivisions; i++) {
        const lineData2 = [];
        for (let j = i; j < summedVectors.length; j += (xDivisions + 1)) {
          lineData2.push(summedVectors[j]);
        }
        data2.datasets.push({
          type: 'line',
          data: lineData2,
          borderColor: 'blue',
          borderWidth: 1,
          fill: false
        });
      }
        
        new Chart(ctx2, {
          type: 'line',
          data: data2,
          options: options2
        });
      
      new Chart(ctx3, {
        type: 'scatter',
        data: data,
        options: options
      });
    }

    const buildButton = document.getElementById('build-button');
    buildButton.addEventListener('click', buildBeam);
  });



  
  
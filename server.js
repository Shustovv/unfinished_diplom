const fs = require('fs');
const http = require('http');
const path = require('path');

http.createServer((req, res) => {
  if (req.url === '/points') {
    const filePath = path.join(__dirname, 'points.json');
    fs.readFile(filePath, 'utf8', (err, data) => {
      if (err) {
        res.statusCode = 500;
        res.end('Error reading file');
      } else {
        res.setHeader('Content-Type', 'application/json');
        res.end(data);
      }
    });
  } else if (req.url === '/') {
    const filePath = path.join(__dirname, 'index.html');
    fs.readFile(filePath, 'utf8', (err, data) => {
      if (err) {
        res.statusCode = 500;
        res.end('Error reading file');
      } else {
        res.setHeader('Content-Type', 'text/html');
        res.end(data);
      }
    });
  } else if (req.url === '/style.css') {
    const filePath = path.join(__dirname, 'style.css');
    fs.readFile(filePath, 'utf8', (err, data) => {
      if (err) {
        res.statusCode = 500;
        res.end('Error reading file');
      } else {
        res.setHeader('Content-Type', 'text/css');
        res.end(data);
      }
    });
  } else if (req.url === '/script.js') {
    const filePath = path.join(__dirname, 'script.js');
    fs.readFile(filePath, 'utf8', (err, data) => {
      if (err) {
        res.statusCode = 500;
        res.end('Error reading file');
      } else {
        res.setHeader('Content-Type', 'application/javascript');
        res.end(data);
      }
    });
  } else {
    res.statusCode = 404;
    res.end('Not found');
  }
}).listen(3000, () => {
  console.log('Server running at http://localhost:3000/');
});

<?php
/**
 * solve_solid.php -- grand-finale solid-viewer API.
 * Pipes a sanitised request to solve_solid.py (the verified solver) and returns
 * its JSON. Value returned == what the verification scripts would compute.
 *
 * Request (POST JSON): { "shape":"cube|ball|torus",
 *   "physics":"heat|wave|liquid|gas", "steps":60, "g":1.0, "c2":0.18 }
 */
header('Content-Type: application/json; charset=utf-8');
header('Access-Control-Allow-Origin: *');
header('Access-Control-Allow-Methods: POST, GET, OPTIONS');
header('Access-Control-Allow-Headers: Content-Type');
if ($_SERVER['REQUEST_METHOD'] === 'OPTIONS') { http_response_code(204); exit; }

$raw = file_get_contents('php://input');
$req = json_decode($raw, true);
if (!is_array($req)) { $req = $_GET; }

$shapes  = ['cube', 'ball', 'torus'];
$physics = ['heat', 'wave', 'liquid', 'gas'];
$shape = isset($req['shape'])   ? (string)$req['shape']   : 'cube';
$phys  = isset($req['physics']) ? (string)$req['physics'] : 'heat';
if (!in_array($shape, $shapes, true) || !in_array($phys, $physics, true)) {
    http_response_code(400);
    echo json_encode(['ok' => false, 'error' => 'invalid shape/physics']);
    exit;
}
$payload = ['shape' => $shape, 'physics' => $phys];
$payload['steps'] = isset($req['steps']) ? max(0, min(2000, intval($req['steps']))) : 60;
$payload['g']     = isset($req['g'])     ? max(0.0, min(50.0, floatval($req['g']))) : 1.0;
$payload['c2']    = isset($req['c2'])    ? max(1e-3, min(1e3, floatval($req['c2']))) : 0.18;

$script = __DIR__ . '/solve_solid.py';
$desc = [0 => ['pipe', 'r'], 1 => ['pipe', 'w'], 2 => ['pipe', 'w']];
$proc = proc_open('python3 ' . escapeshellarg($script), $desc, $pipes, __DIR__);
if (!is_resource($proc)) {
    http_response_code(500);
    echo json_encode(['ok' => false, 'error' => 'failed to start solver']);
    exit;
}
fwrite($pipes[0], json_encode($payload)); fclose($pipes[0]);
$out = stream_get_contents($pipes[1]);
$err = stream_get_contents($pipes[2]);
fclose($pipes[1]); fclose($pipes[2]);
$code = proc_close($proc);
if ($code !== 0 || $out === '') {
    http_response_code(500);
    echo json_encode(['ok' => false, 'error' => 'solver failed', 'detail' => substr($err, 0, 500)]);
    exit;
}
echo $out;

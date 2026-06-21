<?php
/**
 * solve_fluid.php -- V3 fluid-viewer API.
 * Pipes a sanitised request to solve_fluid.py (verified Stage A/B solvers) and
 * returns its time-series JSON. Every frame == a verified solver output.
 *
 * Request (POST JSON): { "scene":"couette|poiseuille|advection",
 *   "U":1.0, "nu":1.0, "G":-1.0, "kappa":0.004, "ux":1.0, "frames":36 }
 */
header('Content-Type: application/json; charset=utf-8');
header('Access-Control-Allow-Origin: *');
header('Access-Control-Allow-Methods: POST, GET, OPTIONS');
header('Access-Control-Allow-Headers: Content-Type');
if ($_SERVER['REQUEST_METHOD'] === 'OPTIONS') { http_response_code(204); exit; }

$raw = file_get_contents('php://input');
$req = json_decode($raw, true);
if (!is_array($req)) { $req = $_GET; }

$scenes = ['couette', 'poiseuille', 'advection'];
$scene = isset($req['scene']) ? (string)$req['scene'] : 'couette';
if (!in_array($scene, $scenes, true)) {
    http_response_code(400);
    echo json_encode(['ok' => false, 'error' => 'invalid scene']);
    exit;
}
$payload = ['scene' => $scene];
$payload['U']      = isset($req['U'])      ? max(-10.0, min(10.0, floatval($req['U'])))      : 1.0;
$payload['nu']     = isset($req['nu'])     ? max(1e-3, min(10.0, floatval($req['nu'])))      : 1.0;
$payload['G']      = isset($req['G'])      ? max(-20.0, min(20.0, floatval($req['G'])))      : -1.0;
$payload['kappa']  = isset($req['kappa'])  ? max(0.0, min(1.0, floatval($req['kappa'])))     : 0.004;
$payload['ux']     = isset($req['ux'])     ? max(-5.0, min(5.0, floatval($req['ux'])))       : 1.0;
$payload['frames'] = isset($req['frames']) ? max(2, min(60, intval($req['frames'])))         : 36;

$script = __DIR__ . '/solve_fluid.py';
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

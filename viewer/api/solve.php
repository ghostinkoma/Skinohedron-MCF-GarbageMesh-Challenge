<?php
/**
 * solve.php -- pressure API.
 *
 * Receives a JSON body (or query params) describing a pressure request, pipes it
 * to the verified Python solver (solve_pressure.py), and returns its JSON. The
 * value returned is exactly what pressure_field_verify.py would compute.
 *
 * Request (POST JSON or GET query):
 *   { "mode":"hydrostatic|atmosphere|acoustic",
 *     "rho":1.0, "g":[0,0,-1.0], "c2":0.15, "steps":120 }
 *
 * Deploy alongside solve_pressure.py. Requires python3 on the server PATH with
 * numpy/scipy and the repo's src/ksf3d importable (the script adds it to sys.path).
 */
header('Content-Type: application/json; charset=utf-8');
header('Access-Control-Allow-Origin: *');
header('Access-Control-Allow-Methods: POST, GET, OPTIONS');
header('Access-Control-Allow-Headers: Content-Type');
if ($_SERVER['REQUEST_METHOD'] === 'OPTIONS') { http_response_code(204); exit; }

// --- gather request ---
$raw = file_get_contents('php://input');
$req = json_decode($raw, true);
if (!is_array($req)) { $req = $_GET; }            // fall back to query params

// --- whitelist + sanitise (never pass arbitrary input straight to the solver) ---
$allowed_modes = ['hydrostatic', 'atmosphere', 'acoustic'];
$mode = isset($req['mode']) ? (string)$req['mode'] : 'hydrostatic';
if (!in_array($mode, $allowed_modes, true)) {
    http_response_code(400);
    echo json_encode(['ok' => false, 'error' => 'invalid mode']);
    exit;
}
$payload = ['mode' => $mode];
$payload['rho'] = isset($req['rho']) ? floatval($req['rho']) : 1.0;
$payload['c2']  = isset($req['c2'])  ? floatval($req['c2'])  : 0.15;
$payload['steps'] = isset($req['steps']) ? intval($req['steps']) : 120;
// gravity vector
if (isset($req['g']) && is_array($req['g']) && count($req['g']) === 3) {
    $payload['g'] = [floatval($req['g'][0]), floatval($req['g'][1]), floatval($req['g'][2])];
} else {
    $payload['g'] = [0.0, 0.0, -1.0];
}
// clamp ranges to keep the solver well-posed
$payload['rho']   = max(-1e3, min(1e3, $payload['rho']));
$payload['c2']    = max(1e-4, min(1e6, $payload['c2']));
$payload['steps'] = max(0, min(2000, $payload['steps']));

// --- invoke the Python solver ---
$python = 'python3';
$script = __DIR__ . '/solve_pressure.py';
$descriptors = [
    0 => ['pipe', 'r'],   // stdin
    1 => ['pipe', 'w'],   // stdout
    2 => ['pipe', 'w'],   // stderr
];
$proc = proc_open(escapeshellarg($python) . ' ' . escapeshellarg($script),
                  $descriptors, $pipes, __DIR__);
if (!is_resource($proc)) {
    http_response_code(500);
    echo json_encode(['ok' => false, 'error' => 'failed to start solver']);
    exit;
}
fwrite($pipes[0], json_encode($payload));
fclose($pipes[0]);
$out = stream_get_contents($pipes[1]);
$err = stream_get_contents($pipes[2]);
fclose($pipes[1]); fclose($pipes[2]);
$code = proc_close($proc);

if ($code !== 0 || $out === '') {
    http_response_code(500);
    echo json_encode(['ok' => false, 'error' => 'solver failed',
                      'detail' => substr($err, 0, 500)]);
    exit;
}
// pass the solver's JSON straight through
echo $out;

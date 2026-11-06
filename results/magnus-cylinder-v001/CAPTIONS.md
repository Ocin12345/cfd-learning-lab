Actual native viscous rotating-cylinder startup, Re_D20, clockwise wall Omega_z=-2 rad/s on a stationary mesh. This full export selects 101 actual states, steps 0 to 500, physical times 0.00 to 10.00 s, at60 ms playback per frame. Near-cylinder view x[-0.35,0.65], y[-0.30,0.30] m uses equal xy aspect, exact128-facet polygon hole and fixed common speed0..0.25 m/s. Saturated native cell counts are retained/labeled. No temporal interpolation, particles, synthetic signed vorticity or forced map uniqueness. Genuine repeated maps remain. Fixed palette affects preview colors only.

Final native gauge-pressure maps share one range from both complete final native fields, with no mean subtraction. Surface Cp uses actual saved/native-boundary pressure at128 native facet centres, outlet0 Pa and q=0.005 Pa; no interpolated or appended wall points. Fluent native face IDs are1-based; OpenFOAM global native face indices are0-based; canonical source face IDs are1-based16257..16384. Any saved OpenFOAM zeroGradient patch without a serialized value uses its exact native owner-pressure semantics, not interpolation. Full-domain uncolored geometry provides context.

Initialization is prescribed rest with actual active inlet/rotating-wall boundary fields, not solved conservation. Unavailable initial ratios stay unavailable; mass/divergence/residual gates and force histories begin at solved step1. Raw native time-zero function outputs remain separate. Pressure/viscous/total force coefficients use each code's actual depth and shared rho/U/D reference. Actual computed lift direction is an outcome, never a gate. Native residual norms differ; cell flux Courant is descriptive. No ideal-flow equality, plateau, periodicity, physical validation or mesh/time independence is claimed.

ANSYS Fluent: primary numerical checks passed
OpenFOAM: primary numerical checks passed

OpenFOAM exact native diagnostic_comparisons:
```json
{
  "primary_prepared_and_original_pilot_inputs_unchanged": true,
  "strict_all_derived_copies_preserved": false,
  "derived_outputs": [
    {
      "path": "0/vorticity",
      "original_sha256": "802235d88b44eb55d7df068170c7ec1b84ef77421baff90a67ea8e62e57add9e",
      "current_sha256": "802235d88b44eb55d7df068170c7ec1b84ef77421baff90a67ea8e62e57add9e",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "0/wallShearStress",
      "original_sha256": "e188989ad36d5472b6dd1c2262d919b992197327acb63f8e7c91ef5ae3c70184",
      "current_sha256": "e188989ad36d5472b6dd1c2262d919b992197327acb63f8e7c91ef5ae3c70184",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0/vorticity",
      "original_sha256": "1c46d8ee960bb83991783b6460017f57262a91a2556b8f142ee26f46ad72c6b2",
      "current_sha256": "1c46d8ee960bb83991783b6460017f57262a91a2556b8f142ee26f46ad72c6b2",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0/wallShearStress",
      "original_sha256": "9f30b51d558848693ae9995eaf044d19a7519092088f76454670021f7032cb32",
      "current_sha256": "9f30b51d558848693ae9995eaf044d19a7519092088f76454670021f7032cb32",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.02/vorticity",
      "original_sha256": "9d52079a65e29f0098fa6cbb7adf19cd1ce245a54e0403694548a2b2caa643d5",
      "current_sha256": "9d52079a65e29f0098fa6cbb7adf19cd1ce245a54e0403694548a2b2caa643d5",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.02/wallShearStress",
      "original_sha256": "cf2488e6a75c6e75d35d5b8cb07bfb331b8a17968622a8a27818b37b7f38b620",
      "current_sha256": "cf2488e6a75c6e75d35d5b8cb07bfb331b8a17968622a8a27818b37b7f38b620",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.04/vorticity",
      "original_sha256": "f904977d93d8900119f0cd13c295459fd93ba8b01e70cc5cc15ffc9ab0474ca5",
      "current_sha256": "f904977d93d8900119f0cd13c295459fd93ba8b01e70cc5cc15ffc9ab0474ca5",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.04/wallShearStress",
      "original_sha256": "2dd0ef204c01f8bf314297f728f98a201ebe79d69b3d8ae6cc08ed41d5450917",
      "current_sha256": "2dd0ef204c01f8bf314297f728f98a201ebe79d69b3d8ae6cc08ed41d5450917",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.06/vorticity",
      "original_sha256": "784d9d07c0afc1ce0b6126b7f5f32a8974b305120ea389ff9db8dd568744b57e",
      "current_sha256": "784d9d07c0afc1ce0b6126b7f5f32a8974b305120ea389ff9db8dd568744b57e",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.06/wallShearStress",
      "original_sha256": "dd85056b69fbdd3d7605b3fc6dcd22719a072fb5745f91fceed5a3130a8272f3",
      "current_sha256": "dd85056b69fbdd3d7605b3fc6dcd22719a072fb5745f91fceed5a3130a8272f3",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.08/vorticity",
      "original_sha256": "b5c31fa347c5a8f454e344ee04362a972dbf51ed9e8459176aab3ad918c072d9",
      "current_sha256": "b5c31fa347c5a8f454e344ee04362a972dbf51ed9e8459176aab3ad918c072d9",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.08/wallShearStress",
      "original_sha256": "6a3db75b09345f65cfc3dccc65300d7cf0e81fc369d643ee60377a33a0b9a6aa",
      "current_sha256": "6a3db75b09345f65cfc3dccc65300d7cf0e81fc369d643ee60377a33a0b9a6aa",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.1/vorticity",
      "original_sha256": "caba14ec849006e69b54eb9ea1d1569eeb41ce86d84e5443d75fc8b27ddcabcf",
      "current_sha256": "7b878440c7eb3bf8ab7535c116cee8acd04a5497ed785ad96c7f375e78a0158e",
      "byte_preserved": false,
      "native_internal_max_abs_difference": 3.799627279477136e-12,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.1/wallShearStress",
      "original_sha256": "b21a0ac5252f20e515b5ac53ac1c85ae348981e62c3a18bcc82cdf6d8cd26b02",
      "current_sha256": "b38b25bb8fa54ec71500cbc5cb73b52071e7b0a63cd3c83f599e020dbdaa3c81",
      "byte_preserved": false,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0/vorticity",
      "original_sha256": "919536a493a80fd38aa0d2141de613db00f2f200cc58920c34eae0021f30ffa2",
      "current_sha256": "919536a493a80fd38aa0d2141de613db00f2f200cc58920c34eae0021f30ffa2",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0/wallShearStress",
      "original_sha256": "9b3d27d47f667be88b64473d1e49e76b48e3dd9aeee76e0fd69338616eb29c32",
      "current_sha256": "9b3d27d47f667be88b64473d1e49e76b48e3dd9aeee76e0fd69338616eb29c32",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.02/vorticity",
      "original_sha256": "bd12e44cefe71415292eb2afe568bfe8ae203a6aabe5500352c3de4a8be82636",
      "current_sha256": "bd12e44cefe71415292eb2afe568bfe8ae203a6aabe5500352c3de4a8be82636",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.02/wallShearStress",
      "original_sha256": "155982e7018e151931a347b1206c076b2bb6ef99d2a0c3241c5cfe93a772b3dc",
      "current_sha256": "155982e7018e151931a347b1206c076b2bb6ef99d2a0c3241c5cfe93a772b3dc",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.04/vorticity",
      "original_sha256": "9044bd334eff3599bd884ab26c4a1acd73e714762e591c5955256fbc971e2179",
      "current_sha256": "9044bd334eff3599bd884ab26c4a1acd73e714762e591c5955256fbc971e2179",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.04/wallShearStress",
      "original_sha256": "05091244fc7d442d24ed88783a0fc39a058a2b2944d0d0d5a981ae8af436a904",
      "current_sha256": "05091244fc7d442d24ed88783a0fc39a058a2b2944d0d0d5a981ae8af436a904",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.06/vorticity",
      "original_sha256": "f385f353690b2294bb42e85932db69dc6face2d87300c714fae48dd06d29733d",
      "current_sha256": "f385f353690b2294bb42e85932db69dc6face2d87300c714fae48dd06d29733d",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.06/wallShearStress",
      "original_sha256": "8fa091d3a2c659d02829fb84ee54de650cba34cafa056c7a1c5483de9e257564",
      "current_sha256": "8fa091d3a2c659d02829fb84ee54de650cba34cafa056c7a1c5483de9e257564",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.08/vorticity",
      "original_sha256": "88a2ba56debd1bcdcd361a6a2cb0ff26944dc3383aa7c7ef8e642603d9444afe",
      "current_sha256": "88a2ba56debd1bcdcd361a6a2cb0ff26944dc3383aa7c7ef8e642603d9444afe",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.08/wallShearStress",
      "original_sha256": "26585fb845be212c3186dd706fa0db8be171cf606ba30756873ef7c5e01057e5",
      "current_sha256": "26585fb845be212c3186dd706fa0db8be171cf606ba30756873ef7c5e01057e5",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.1/vorticity",
      "original_sha256": "8ba009879e6f4268d6b79d8bdcd4280a0fef045168c568401a32ee3712ab0498",
      "current_sha256": "e103ba15587a463ed8e46a22afe661a2ea32bc025f73b32dd86f7ba9c254f266",
      "byte_preserved": false,
      "native_internal_max_abs_difference": 4.299671729768306e-12,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.1/wallShearStress",
      "original_sha256": "328b5d85fd63c01840ba4aba4dee5d3ed711facbb0ab9e5348cef5e398d2f68e",
      "current_sha256": "76158f07fff169cc10d5ec3b33ea1988ba178e8a0594d828090f16ccb280b0ce",
      "byte_preserved": false,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    }
  ],
  "original_pilot_diagnostics": {
    "primary_initial_inputs_unchanged": true,
    "strict_all_derived_copies_preserved": true,
    "derived_outputs": [
      {
        "path": "0/vorticity",
        "original_sha256": "802235d88b44eb55d7df068170c7ec1b84ef77421baff90a67ea8e62e57add9e",
        "current_sha256": "802235d88b44eb55d7df068170c7ec1b84ef77421baff90a67ea8e62e57add9e",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      },
      {
        "path": "0/wallShearStress",
        "original_sha256": "e188989ad36d5472b6dd1c2262d919b992197327acb63f8e7c91ef5ae3c70184",
        "current_sha256": "e188989ad36d5472b6dd1c2262d919b992197327acb63f8e7c91ef5ae3c70184",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      },
      {
        "path": "processor0/0/vorticity",
        "original_sha256": "1c46d8ee960bb83991783b6460017f57262a91a2556b8f142ee26f46ad72c6b2",
        "current_sha256": "1c46d8ee960bb83991783b6460017f57262a91a2556b8f142ee26f46ad72c6b2",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      },
      {
        "path": "processor0/0/wallShearStress",
        "original_sha256": "9f30b51d558848693ae9995eaf044d19a7519092088f76454670021f7032cb32",
        "current_sha256": "9f30b51d558848693ae9995eaf044d19a7519092088f76454670021f7032cb32",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      },
      {
        "path": "processor1/0/vorticity",
        "original_sha256": "919536a493a80fd38aa0d2141de613db00f2f200cc58920c34eae0021f30ffa2",
        "current_sha256": "919536a493a80fd38aa0d2141de613db00f2f200cc58920c34eae0021f30ffa2",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      },
      {
        "path": "processor1/0/wallShearStress",
        "original_sha256": "9b3d27d47f667be88b64473d1e49e76b48e3dd9aeee76e0fd69338616eb29c32",
        "current_sha256": "9b3d27d47f667be88b64473d1e49e76b48e3dd9aeee76e0fd69338616eb29c32",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      }
    ]
  },
  "original_pilot0_to5_authoritative": true
}
```

OpenFOAM exact native execution_completion:
```json
{
  "actual_solver_completed": true,
  "inherited_actual_steps": 5,
  "additional_actual_steps": 495,
  "total_actual_steps": 500,
  "first_new_time_s": 0.12,
  "first_new_index": 6,
  "native_rank_states": 501,
  "rank_time_records": 1000,
  "raw_wrapper_status": {
    "stage": "finished_production_execution",
    "solver_started": true,
    "mpi_ranks": 2,
    "inherited_actual_steps": 5,
    "additional_actual_steps": 495,
    "total_actual_steps": 500,
    "physical_time_s": 10,
    "deltaT_s": 0.02,
    "native_rank_states": 501,
    "outer_correctors": 30,
    "serial_reconstruction_requested": false,
    "production_admitted": true,
    "numerical_analysis_pending": true
  },
  "serial_reconstruction_requested": false
}
```

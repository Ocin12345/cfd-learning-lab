Actual native offset tandem-square velocity fields, two stationary 0.1 m squares, Re_D100. This full export selects 101 actual saved states, steps 0 to 800, physical times 0.00 to 40.00 s, at 60 ms per playback frame. The exact 9472-cell, 1.48 m2 fluid mesh has two solid holes. Maps use equal xy aspect and a fixed common 0 to 0.2 m/s speed scale; any saturation is labeled and actual counts retained. No temporal interpolation, particles, synthesized vorticity or forced map uniqueness. Genuine settled/quantized repeats are preserved. Fixed palette quantization affects preview colors only; static maps and CSV retain native values.

Initialization is prescribed cell rest with actual native boundary flux, not a solved conservation state. Undefined initial relative boundary ratio remains unavailable, never replaced with zero. Mass/divergence histories and residual targets begin at solved step1; force histories include solved steps only. Raw OpenFOAM startup rows at time zero remain separate from solved histories. Per-body pressure, viscous and total coefficients use each code's native depth and the same reference D/U/rho; different-depth raw forces are not compared. Native residual norms remain separate. Courant number is descriptive. Profiles use actual columns ix52/84/112: one between the squares and two downstream, at native cell centres without appended wall values. Native gauge pressure is exported in Pa without mean subtraction.

No steady-force plateau, zero mean lift, periodicity, benchmark agreement, physical validation, mesh independence or time-step independence is claimed. Full native field CSV remains external; compact final maps and all solved-step history rows are provided for the record.

ANSYS Fluent: primary numerical checks passed
OpenFOAM: primary numerical checks passed

OpenFOAM exact native diagnostic_comparisons (separate from primary numerical checks):
```json
{
  "primary_prepared_and_original_pilot_inputs_unchanged": true,
  "strict_all_derived_copies_preserved": false,
  "derived_outputs": [
    {
      "path": "0/vorticity",
      "original_sha256": "ea7654a69c429f836e10a2aa2f51d3b9e32b7f309de8370f8b0bdfeca272caa8",
      "current_sha256": "ea7654a69c429f836e10a2aa2f51d3b9e32b7f309de8370f8b0bdfeca272caa8",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "0/wallShearStress",
      "original_sha256": "05067024d702797ce9747280750adb3cd32850e0cb3456fc0700ca76b17669f6",
      "current_sha256": "05067024d702797ce9747280750adb3cd32850e0cb3456fc0700ca76b17669f6",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0/vorticity",
      "original_sha256": "41f8066630e99d7e30da261d5bf62b87ce4821be836bd96ee58ec550bb8f053f",
      "current_sha256": "41f8066630e99d7e30da261d5bf62b87ce4821be836bd96ee58ec550bb8f053f",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0/wallShearStress",
      "original_sha256": "ea6dc227047d8e5abdc1787c5e2104c2c19c37c26378697a045f7cfa6420f353",
      "current_sha256": "ea6dc227047d8e5abdc1787c5e2104c2c19c37c26378697a045f7cfa6420f353",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.05/vorticity",
      "original_sha256": "a1aba4da19dc158d26cdc9d014be0373c5fd2c477ae25c0146527ffbc06434ea",
      "current_sha256": "a1aba4da19dc158d26cdc9d014be0373c5fd2c477ae25c0146527ffbc06434ea",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.05/wallShearStress",
      "original_sha256": "44cd4437ad0fb790021df48f5c1c3773db0988f8723462a60ec7ca64a395c2da",
      "current_sha256": "44cd4437ad0fb790021df48f5c1c3773db0988f8723462a60ec7ca64a395c2da",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.1/vorticity",
      "original_sha256": "e558fb229df75e1914475c1b23b69b0a2cfa84d059ccf638ffd760d0bb228da4",
      "current_sha256": "e558fb229df75e1914475c1b23b69b0a2cfa84d059ccf638ffd760d0bb228da4",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.1/wallShearStress",
      "original_sha256": "abb1c5d07812483726c4ff5de1bfadd13d463e586f785f18496ee9ff5787618b",
      "current_sha256": "abb1c5d07812483726c4ff5de1bfadd13d463e586f785f18496ee9ff5787618b",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.15/vorticity",
      "original_sha256": "4b5be806a6a5ac8d20dfbd4cb67966afff99a8b56c31bcc57bed20988697c73a",
      "current_sha256": "4b5be806a6a5ac8d20dfbd4cb67966afff99a8b56c31bcc57bed20988697c73a",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.15/wallShearStress",
      "original_sha256": "9f85957f2e707b171c837dc4e2deb1bfdadc1513d3be323431c1c280a61580bf",
      "current_sha256": "9f85957f2e707b171c837dc4e2deb1bfdadc1513d3be323431c1c280a61580bf",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.2/vorticity",
      "original_sha256": "a00857abd2dab1cf48ae2027ec46c722f749a9d46d16b397118dfd33147b846d",
      "current_sha256": "a00857abd2dab1cf48ae2027ec46c722f749a9d46d16b397118dfd33147b846d",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.2/wallShearStress",
      "original_sha256": "460e9099cb189f721a5bf305a3c95f41bd80287904a2931dd146f25415827c4d",
      "current_sha256": "460e9099cb189f721a5bf305a3c95f41bd80287904a2931dd146f25415827c4d",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.25/vorticity",
      "original_sha256": "f115134a3a8df4218f561f35a803f25faaf2e69cbeca351590fef62ef61a859d",
      "current_sha256": "2f4202de60147ec0dfc1fed9f0b666d9f5c85114184130e57b1a789b3c52837b",
      "byte_preserved": false,
      "native_internal_max_abs_difference": 1.000088900582341e-12,
      "both_fields_finite": true
    },
    {
      "path": "processor0/0.25/wallShearStress",
      "original_sha256": "df2f558764623a7d7ace5387fbe3e98b228db6c26894749fff68782b97a15f68",
      "current_sha256": "a5850271d94800856cb40b9ccf8949627c135617211c099d8edd1295ec27618d",
      "byte_preserved": false,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0/vorticity",
      "original_sha256": "3decd0a949bacd943dbdd6d5a2295c1428736d5efb10c0ba2ae0f9970fe18102",
      "current_sha256": "3decd0a949bacd943dbdd6d5a2295c1428736d5efb10c0ba2ae0f9970fe18102",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0/wallShearStress",
      "original_sha256": "e555f6492a8f1f0548a546b8f6ccc7a0eb6449f623d87f77f5d095674899c1a0",
      "current_sha256": "e555f6492a8f1f0548a546b8f6ccc7a0eb6449f623d87f77f5d095674899c1a0",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.05/vorticity",
      "original_sha256": "f41a0d65954bc135118512755130595393d3adb3e89810711605975708a0825a",
      "current_sha256": "f41a0d65954bc135118512755130595393d3adb3e89810711605975708a0825a",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.05/wallShearStress",
      "original_sha256": "abca582e992cb2c31c6223aec42ca48e4845df98a3b5928d1102d9fce5b475c8",
      "current_sha256": "abca582e992cb2c31c6223aec42ca48e4845df98a3b5928d1102d9fce5b475c8",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.1/vorticity",
      "original_sha256": "51487116800482866202dfd6fb87b512bd2f5b27aa787def2d9dccf49cebdff3",
      "current_sha256": "51487116800482866202dfd6fb87b512bd2f5b27aa787def2d9dccf49cebdff3",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.1/wallShearStress",
      "original_sha256": "14d135062aee4914e62732aed1f97f6edb312dea085783c8615f5019237e69e6",
      "current_sha256": "14d135062aee4914e62732aed1f97f6edb312dea085783c8615f5019237e69e6",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.15/vorticity",
      "original_sha256": "9cca16b100f5c2461daa6bb11d05da7c5c45f4e1e1b15e439281738312b85100",
      "current_sha256": "9cca16b100f5c2461daa6bb11d05da7c5c45f4e1e1b15e439281738312b85100",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.15/wallShearStress",
      "original_sha256": "cff080c4533f887eb1915a62c3dd74536b50818021d2267d5ce45dcc9c221a10",
      "current_sha256": "cff080c4533f887eb1915a62c3dd74536b50818021d2267d5ce45dcc9c221a10",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.2/vorticity",
      "original_sha256": "e69f943b82022eaa04630490e4015334065d6281e720aeadbbf827b2b71ece58",
      "current_sha256": "e69f943b82022eaa04630490e4015334065d6281e720aeadbbf827b2b71ece58",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.2/wallShearStress",
      "original_sha256": "4c8e0c3c332568feb245655a969aacf6d7ef02ab49fa7d3115109a81dec1fe04",
      "current_sha256": "4c8e0c3c332568feb245655a969aacf6d7ef02ab49fa7d3115109a81dec1fe04",
      "byte_preserved": true,
      "native_internal_max_abs_difference": 0.0,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.25/vorticity",
      "original_sha256": "5dd4ecfef59a8a66ce3fbade66c5f7d231c17effe9107773e24f8be6aaf41aaf",
      "current_sha256": "f35d9e04f8c6b15360c85cc3339a71946c723ced6b56555ba0797b92fb291990",
      "byte_preserved": false,
      "native_internal_max_abs_difference": 3.960482115799921e-13,
      "both_fields_finite": true
    },
    {
      "path": "processor1/0.25/wallShearStress",
      "original_sha256": "1877e7d6c2d76492abf2e6c1fc8166c58491aca3bf7488a89fcd5682b0f97e8b",
      "current_sha256": "1877e7d6c2d76492abf2e6c1fc8166c58491aca3bf7488a89fcd5682b0f97e8b",
      "byte_preserved": true,
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
        "original_sha256": "ea7654a69c429f836e10a2aa2f51d3b9e32b7f309de8370f8b0bdfeca272caa8",
        "current_sha256": "ea7654a69c429f836e10a2aa2f51d3b9e32b7f309de8370f8b0bdfeca272caa8",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      },
      {
        "path": "0/wallShearStress",
        "original_sha256": "05067024d702797ce9747280750adb3cd32850e0cb3456fc0700ca76b17669f6",
        "current_sha256": "05067024d702797ce9747280750adb3cd32850e0cb3456fc0700ca76b17669f6",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      },
      {
        "path": "processor0/0/vorticity",
        "original_sha256": "41f8066630e99d7e30da261d5bf62b87ce4821be836bd96ee58ec550bb8f053f",
        "current_sha256": "41f8066630e99d7e30da261d5bf62b87ce4821be836bd96ee58ec550bb8f053f",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      },
      {
        "path": "processor0/0/wallShearStress",
        "original_sha256": "ea6dc227047d8e5abdc1787c5e2104c2c19c37c26378697a045f7cfa6420f353",
        "current_sha256": "ea6dc227047d8e5abdc1787c5e2104c2c19c37c26378697a045f7cfa6420f353",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      },
      {
        "path": "processor1/0/vorticity",
        "original_sha256": "3decd0a949bacd943dbdd6d5a2295c1428736d5efb10c0ba2ae0f9970fe18102",
        "current_sha256": "3decd0a949bacd943dbdd6d5a2295c1428736d5efb10c0ba2ae0f9970fe18102",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      },
      {
        "path": "processor1/0/wallShearStress",
        "original_sha256": "e555f6492a8f1f0548a546b8f6ccc7a0eb6449f623d87f77f5d095674899c1a0",
        "current_sha256": "e555f6492a8f1f0548a546b8f6ccc7a0eb6449f623d87f77f5d095674899c1a0",
        "byte_preserved": true,
        "native_internal_max_abs_difference": 0.0,
        "both_fields_finite": true
      }
    ]
  },
  "original_pilot0_to5_authoritative": true
}
```

OpenFOAM exact native execution_completion (separate from primary numerical checks):
```json
{
  "actual_solver_completed": true,
  "inherited_actual_steps": 5,
  "additional_actual_steps": 795,
  "total_actual_steps": 800,
  "first_new_time_s": 0.3,
  "first_new_index": 6,
  "native_rank_states": 801,
  "rank_time_records": 1600,
  "raw_wrapper_status": {
    "stage": "finished_production_execution",
    "solver_started": true,
    "mpi_ranks": 2,
    "inherited_actual_steps": 5,
    "additional_actual_steps": 795,
    "total_actual_steps": 800,
    "physical_time_s": 40,
    "deltaT_s": 0.05,
    "native_rank_states": 801,
    "outer_correctors": 30,
    "serial_reconstruction_requested": false,
    "production_admitted": true,
    "numerical_analysis_pending": true
  },
  "serial_reconstruction_requested": false
}
```

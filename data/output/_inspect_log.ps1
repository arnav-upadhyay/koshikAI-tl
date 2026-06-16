$csv = Import-Csv 'e:/TxTl_v2_main/TXTLv2/data/output/spindler_rplL_elongation_log.csv'
$csv | Group-Object timestamp, N | ForEach-Object {
  $f = $_.Group[0]
  [PSCustomObject]@{
    ts    = $f.timestamp
    N     = [int]$f.N
    Stau  = $f.rplL_total_time_s
    mv    = $f.rplL_mean_v_cs
    T3m   = $f.T3_mean
    mbg   = $f.m_bg_each
    VIN   = $f.V_IN_70SIC_scale
    ARS   = $f.ARS_scale
    tRNA  = $f.tRNA_scale
    EF    = $f.EF_scale
    Rib   = $f.Rib_scale
    atol  = $f.ode_abstol
    wall  = $f.wall_time_s
    rows  = $_.Count
  }
} | Sort-Object N, ts | Format-Table -AutoSize

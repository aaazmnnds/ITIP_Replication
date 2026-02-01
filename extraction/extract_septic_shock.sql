-- extract_septic_shock.sql
-- Reference extraction script for MIMIC-III Septic Shock Cohort (Sepsis-3)
-- Based on the official MIMIC Code Repository logic (MIT-LCP/mimic-code)

-- 1. Identify Suspected Infection
-- Patients with antibiotic administration within +/- 72 hours of a microbiology culture.
WITH suspected_infection AS (
  SELECT 
    subject_id, hadm_id, stay_id,
    prob_infection_time
  FROM `physionet-data.mimic_derived.suspected_infection`
),

-- 2. Calculate SOFA Score
-- Extract max SOFA score within the infection window.
sofa_scores AS (
  SELECT 
    stay_id,
    hr,
    respiration_24hours as sofa_respiration,
    coagulation_24hours as sofa_coagulation,
    liver_24hours as sofa_liver,
    cardiovascular_24hours as sofa_cardiovascular,
    cns_24hours as sofa_cns,
    renal_24hours as sofa_renal,
    sofa_24hours as sofa_score
  FROM `physionet-data.mimic_derived.sofa`
),

-- 3. Identify Sepsis-3 Patients
-- Infection + SOFA >= 2
sepsis3 AS (
  SELECT 
    si.subject_id, si.hadm_id, si.stay_id,
    si.prob_infection_time,
    s.sofa_score
  FROM suspected_infection si
  JOIN sofa_scores s ON si.stay_id = s.stay_id
  WHERE s.sofa_score >= 2
),

-- 4. Identify Septic Shock
-- Sepsis-3 + Vasopressor Use + Lactate > 2 mmol/L (despite fluid resuscitation)
septic_shock AS (
  SELECT 
    s3.subject_id, s3.hadm_id, s3.stay_id,
    s3.prob_infection_time,
    v.vasopressor_time,
    l.valuenum as lactate_level
  FROM sepsis3 s3
  LEFT JOIN `physionet-data.mimic_derived.vasopressor_duration` v
    ON s3.stay_id = v.stay_id
  LEFT JOIN `physionet-data.mimic_hosp.labevents` l
    ON s3.hadm_id = l.hadm_id
    AND l.itemid IN (50813) -- Lactate
  WHERE v.vasopressor_time IS NOT NULL
    AND l.valuenum > 2.0
),

-- 5. Extract Feature Set (Vitals, Labs, Demographics)
feature_extraction AS (
  SELECT
    ss.subject_id, ss.hadm_id, ss.stay_id,
    -- Demographics
    p.gender, p.anchor_age as age,
    
    -- Vitals (Mean within first 24h of shock onset)
    AVG(v.heart_rate) as heart_rate_mean,
    AVG(v.sbp) as sysbp_mean,
    AVG(v.dbp) as diasbp_mean,
    AVG(v.mbp) as meanbp_mean,
    AVG(v.resp_rate) as resprate_mean,
    AVG(v.temperature) as tempc_mean,
    AVG(v.spo2) as spo2_mean,
    
    -- Labs (First value within 24h)
    -- Grouping not shown for brevity, implies joining with labevents
    MAX(CASE WHEN l.itemid = 50813 THEN l.valuenum END) as lactate_max,
    MIN(CASE WHEN l.itemid = 50912 THEN l.valuenum END) as creatinine_min,
    MIN(CASE WHEN l.itemid = 51006 THEN l.valuenum END) as bun_min,
    MIN(CASE WHEN l.itemid = 50931 THEN l.valuenum END) as glucose_min,
    MIN(CASE WHEN l.itemid = 51221 THEN l.valuenum END) as hematocrit_min,
    MIN(CASE WHEN l.itemid = 51301 THEN l.valuenum END) as wbc_min,
    MIN(CASE WHEN l.itemid = 50821 THEN l.valuenum END) as po2_min,
    MIN(CASE WHEN l.itemid = 50818 THEN l.valuenum END) as pco2_min,
    MIN(CASE WHEN l.itemid = 50820 THEN l.valuenum END) as ph_min,
    
    -- Outcome
    a.hospital_expire_flag as mortality
    
  FROM septic_shock ss
  JOIN `physionet-data.mimic_core.patients` p ON ss.subject_id = p.subject_id
  JOIN `physionet-data.mimic_core.admissions` a ON ss.hadm_id = a.hadm_id
  LEFT JOIN `physionet-data.mimic_icu.chartevents` v ON ss.stay_id = v.stay_id
  LEFT JOIN `physionet-data.mimic_hosp.labevents` l ON ss.hadm_id = l.hadm_id
  
  GROUP BY ss.subject_id, ss.hadm_id, ss.stay_id, p.gender, p.anchor_age, a.hospital_expire_flag
)

SELECT * FROM feature_extraction;

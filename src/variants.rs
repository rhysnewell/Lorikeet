//use std::cmp::Ordering;
//use std::fmt::Debug;
//use std::ops::{Deref, Range};
//
//
//use bio::stats::LogProb;
//use itertools::Itertools;
//use ordered_float::NotNan;
//use strum_macros::{EnumIter, EnumString, IntoStaticStr};
//
//// Variant data type adapted from varlociraptor
//// Only taken because compiling varlociraptor takes a very long time
//// Copyright 2016-2019 Johannes Köster, David Lähnemann.
//#[derive(Debug, Clone, IntoStaticStr)]
//pub enum VariantType {
//    Insertion(Option<Range<u32>>),
//    Deletion(Option<Range<u32>>),
//    SNV,
//    MNV,
//    None, // site with no suggested alternative allele
//}
//
//impl From<&str> for VariantType {
//    fn from(string: &str) -> VariantType {
//        match string {
//            "INS" => VariantType::Insertion(None),
//            "DEL" => VariantType::Deletion(None),
//            "SNV" => VariantType::SNV,
//            "REF" => VariantType::None,
//            _ => panic!("bug: given string does not describe a valid variant type"),
//        }
//    }
//}
//
//#[derive(Clone, Debug)]
//pub enum Variant {
//    Deletion(u32),
//    Insertion(Vec<u8>),
//    SNV(u8),
//    MNV(Vec<u8>),
//    None,
//}
//
//impl Variant {
//    pub fn has_fragment_evidence(&self) -> bool {
//        match self {
//            &Variant::Deletion(_) => true,
//            &Variant::Insertion(_) => true,
//            &Variant::SNV(_) => false,
//            &Variant::MNV(_) => false,
//            &Variant::None => false,
//        }
//    }
//
//    pub fn is_single_base(&self) -> bool {
//        match self {
//            &Variant::SNV(_) | &Variant::None => true,
//            _ => false,
//        }
//    }
//
//    pub fn is_snv(&self) -> bool {
//        match self {
//            &Variant::SNV(_) => true,
//            _ => false,
//        }
//    }
//
//    pub fn is_indel(&self) -> bool {
//        match self {
//            &Variant::Deletion(_) => true,
//            &Variant::Insertion(_) => true,
//            &Variant::SNV(_) => false,
//            &Variant::MNV(_) => false,
//            &Variant::None => false,
//        }
//    }
//
//    pub fn is_type(&self, vartype: &VariantType) -> bool {
//        match (self, vartype) {
//            (&Variant::Deletion(l), &VariantType::Deletion(Some(ref range))) => {
//                l >= range.start && l < range.end
//            }
//            (&Variant::Insertion(_), &VariantType::Insertion(Some(ref range))) => {
//                self.len() >= range.start && self.len() < range.end
//            }
//            (&Variant::Deletion(_), &VariantType::Deletion(None)) => true,
//            (&Variant::Insertion(_), &VariantType::Insertion(None)) => true,
//            (&Variant::SNV(_), &VariantType::SNV) => true,
//            (&Variant::MNV(_), &VariantType::MNV) => true,
//            (&Variant::None, &VariantType::None) => true,
//            _ => false,
//        }
//    }
//
//    pub fn end(&self, start: u32) -> u32 {
//        match self {
//            &Variant::Deletion(length) => start + length,
//            &Variant::Insertion(_) => start + 1, // end of insertion is the next regular base
//            &Variant::SNV(_) | &Variant::None => start,
//            &Variant::MNV(ref alt) => start + alt.len() as u32,
//        }
//    }
//
//    pub fn centerpoint(&self, start: u32) -> u32 {
//        match self {
//            &Variant::Deletion(length) => start + length / 2,
//            &Variant::Insertion(_) => start, // end of insertion is the next regular base
//            &Variant::SNV(_) | &Variant::None => start,
//            &Variant::MNV(ref alt) => start + alt.len() as u32 / 2,
//        }
//    }
//
//    pub fn len(&self) -> u32 {
//        match self {
//            &Variant::Deletion(l) => l,
//            &Variant::Insertion(ref s) => s.len() as u32,
//            &Variant::SNV(_) => 1,
//            &Variant::MNV(ref alt) => alt.len() as u32,
//            &Variant::None => 1,
//        }
//    }
//}

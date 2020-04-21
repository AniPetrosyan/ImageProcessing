var Student = {
  // please fill in your name and NetID
  //  your NetID is the part of your email before @princeton.edu
  name: "Ani_Petrosian",
  tumoID: "ani.petrosyan.y",
};

Student.updateHTML = function () {
  var studentInfo = this.name + " &lt;" + this.tumoID + "&gt;";
  document.getElementById("student").innerHTML = studentInfo;
};
